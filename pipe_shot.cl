#include "common.cl"

struct pipe_specific {
    int seg_no;
    struct id_pipe *pipe_id;
};

struct id_pipe {
    struct lin_pipe *linear;
    struct bend_pipe *curved;
    int pipe_is_bend;
};


struct lin_pipe {
    struct id_pipe* l;
    int pipe_is_bend;
    double pipe_V[3];              /* start point for pipe section */
    double pipe_H[3];              /* unit vector in direction of pipe section */
    double pipe_ribase, pipe_ritop;        /* base and top inner radii */
    double pipe_ribase_sq, pipe_ritop_sq;  /* inner radii squared */
    double pipe_ridiff_sq, pipe_ridiff;    /* difference between top and base inner radii */
    double pipe_rodiff_sq, pipe_rodiff;    /* difference between top and base outer radii */
    double pipe_robase, pipe_rotop;        /* base and top outer radii */
    double pipe_robase_sq, pipe_rotop_sq;  /* outer radii squared */
    double pipe_len;               /* length of pipe segment */
    double pipe_SoR[16];             /* Scale and rotate */
    double pipe_invRoS[16];              /* inverse rotation and scale */
    double pipe_min[3];
    double pipe_max[3];
};


struct bend_pipe {
    struct id_pipe* l;
    int pipe_is_bend;
    double bend_radius;        /* distance from bend_v to center of pipe */
    double bend_or;            /* outer radius */
    double bend_ir;            /* inner radius */
    double bend_invR[16];            /* inverse rotation matrix */
    double bend_SoR[16];         /* Scale and rotate */
    double bend_V[3];         /* Center of bend */
    double bend_start[3];         /* Start of bend */
    double bend_end[3];           /* End of bend */
    double bend_alpha_i;       /* ratio of inner radius to bend radius */
    double bend_alpha_o;       /* ratio of outer radius to bend radius */
    double bend_angle;         /* Angle that bend goes through */
    double bend_ra[3];         /* unit vector in plane of bend (points toward start from bend_V) */
    double bend_rb[3];         /* unit vector in plane of bend (normal to bend_ra) */
    double bend_endNorm[3];        /* unit vector normal to end plane */
    double bend_startNorm[3];      /* unit vector normal to start plane */
    double bend_N[3];          /* unit vector normal to plane of bend */
    double bend_bound_center[3];      /* center of bounding sphere */
    double bend_bound_radius_sq;   /* square of bounding sphere radius */
};



#define PIPE_CONNECTING_ARCS 4 /* number of connecting arcs to draw between points */
#define PIPE_CIRCLE_SEGS 16    /* number of segments used to plot a circle */

#define PIPE_LINEAR_OUTER_BODY  1
#define PIPE_LINEAR_INNER_BODY  2
#define PIPE_LINEAR_TOP     3
#define PIPE_LINEAR_BASE    4
#define PIPE_BEND_OUTER_BODY    5
#define PIPE_BEND_INNER_BODY    6
#define PIPE_BEND_BASE      7
#define PIPE_BEND_TOP       8
#define PIPE_RADIUS_CHANGE      9

#define RT_PIPE_MAXHITS 128



/**
 * Check for hits on surfaces created by discontinuous radius changes
 * from one pipe segment to the next. Can only happen when one segment
 * is a bend, because linear segments handle different radii at each
 * end. Bend segments must have constant radii .  These surfaces are
 * normal to the flow of the pipe.
 */
 void discont_radius_shot(double3 r_pt, double3 r_dir, double3 center, double3 norm, double or1_sq, double ir1_sq, 
    double or2_sq, double ir2_sq, struct hit *hits, int *hit_count, int seg_no)
{
    double dist_to_plane;
    double norm_dist;
    double slant_factor;
    double t_tmp;
    double3 hit_pt;
    double radius_sq;

    /* calculate intersection with plane at center (with normal "norm") */
    dist_to_plane = dot(norm, center);
    norm_dist = dist_to_plane - dot(norm, r_pt);
    slant_factor = dot(norm, r_dir);
    if (!ZERO(slant_factor)) {
    double3 to_center;
    struct hit *hitp;

    t_tmp = norm_dist / slant_factor;
    hit_pt = r_pt + t_tmp * r_dir;
    to_center = center - hit_pt;
    radius_sq = dot(to_center, to_center);

    /* where the radius ranges overlap, there is no hit */
    if (radius_sq <= or1_sq && radius_sq >= ir1_sq && radius_sq <= or2_sq && radius_sq >= ir2_sq)
        return;
    
    /* if we are within one of the radius ranges, we have a hit */
    if ((radius_sq <= or2_sq && radius_sq >= ir2_sq) || (radius_sq <= or1_sq && radius_sq >= ir1_sq)) 
    {
        hitp = &hits[*hit_count];
        hitp->hit_dist = t_tmp;
        hitp->hit_surfno = seg_no * 10 + PIPE_RADIUS_CHANGE;

        /* within first range, use norm, otherwise reverse */
        if (radius_sq <= or1_sq && radius_sq >= ir1_sq)
            hitp->hit_normal = vload3(0, norm); 
        else 
            hitp->hit_normal = -vload3(0, norm);
        
        if ((*hit_count)++ >= RT_PIPE_MAXHITS)
            return;  
    }
    }
}



 void bend_pipe_shot(double3 r_pt, double3 r_dir, struct bend_pipe *bp, struct hit *hits, int *hit_count, int seg_no)
{
    double3 dprime;      /* D' */
    double3 pprime;      /* P' */
    double3 work;        /* temporary vector */
    double C[5];        /* The final equation */
    bn_complex_t val[4];    /* The complex roots */
    int j;

    int root_count = 0;
    double A[3], Asqr[5];
    double X2_Y2[3];        /* X**2 + Y**2 */
    double3 cor_pprime;      /* new ray origin */
    double cor_proj;
    double or_sq;      /* outside radius squared */
    double ir_sq;      /* inside radius squared */
    double or2_sq;     /* outside radius squared (from adjacent seg) */
    double ir2_sq;     /* inside radius squared (from adjacent seg) */
    int parallel;       /* set to one when ray is parallel to plane of bend */
    double dist = 0;       /* distance between ray and plane of bend */
    double tmp;
    struct id_pipe *prev;
    struct id_pipe *next;

    int ctr=0;
    while(ctr<1)
    {
        or_sq = bp->bend_or * bp->bend_or;
        ir_sq = bp->bend_ir * bp->bend_ir;

        tmp = dot(r_dir, vload3(0, bp->bend_N));
        if (NEAR_ZERO(tmp, 0.0000005)) 
        {
        /* ray is parallel to plane of bend */
        parallel = 1;
        dist = fabs(dot(r_pt, vload3(0, bp->bend_N)) - dot(vload3(0, bp->bend_V), vload3(0, bp->bend_N)));

        if (dist > bp->bend_or) 
            break; /* ray is more than outer radius away from plane of bend */
        
        } 
        else 
            parallel = 0;
        

        /* Convert vector into the space of the unit torus */
        dprime = MAT4X3VEC(bp->bend_SoR, r_dir);
        dprime = normalize(dprime);

        work = r_pt - vload3(0, bp->bend_V);
        pprime = MAT4X3VEC(bp->bend_SoR, work);

        /* normalize distance from torus.  substitute corrected pprime
         * which contains a translation along ray direction to closest
         * approach to vertex of torus.  Translating ray origin along
         * direction of ray to closest pt. to origin of solid's coordinate
         * system, new ray origin is 'cor_pprime'.
         */
        cor_proj = dot(pprime, dprime);
        cor_pprime = dprime * cor_proj;
        cor_pprime = pprime - cor_pprime;

        /*
         * Given a line and a ratio, alpha, finds the equation of the unit
         * torus in terms of the variable 't'.
         *
         * The equation for the torus is:
         *
         * [ X**2 + Y**2 + Z**2 + (1 - alpha**2) ]**2 - 4*(X**2 + Y**2) = 0
         *
         * First, find X, Y, and Z in terms of 't' for this line, then
         * substitute them into the equation above.
         *
         * Wx = Dx*t + Px
         *
         * Wx**2 = Dx**2 * t**2  +  2 * Dx * Px  +  Px**2
         *         [0]                 [1]           [2]    dgr=2
         */

        X2_Y2[0] = dprime.x * dprime.x + dprime.y * dprime.y;
        X2_Y2[1] = 2.0 * (dprime.x * cor_pprime.x + dprime.y * cor_pprime.y);
        X2_Y2[2] = cor_pprime.x * cor_pprime.x + cor_pprime.y * cor_pprime.y;

        /* A = X2_Y2 + Z2 */
        A[0] = X2_Y2[0] + dprime.z * dprime.z;
        A[1] = X2_Y2[1] + 2.0 * dprime.z * cor_pprime.z;
        A[2] = X2_Y2[2] + cor_pprime.z * cor_pprime.z + 1.0 - bp->bend_alpha_o * bp->bend_alpha_o;

        /* Inline expansion of (void) bn_poly_mul(&Asqr, &A, &A) */
        /* Both polys have degree two */
        Asqr[0] = A[0] * A[0];
        Asqr[1] = A[0] * A[1] + A[1] * A[0];
        Asqr[2] = A[0] * A[2] + A[1] * A[1] + A[2] * A[0];
        Asqr[3] = A[1] * A[2] + A[2] * A[1];
        Asqr[4] = A[2] * A[2];

        /* Inline expansion of bn_poly_scale(&X2_Y2, 4.0) and
         * bn_poly_sub(&C, &Asqr, &X2_Y2).
         */
        C[0] = Asqr[0];
        C[1] = Asqr[1];
        C[2] = Asqr[2] - X2_Y2[0] * 4.0;
        C[3] = Asqr[3] - X2_Y2[1] * 4.0;
        C[4] = Asqr[4] - X2_Y2[2] * 4.0;

        /* It is known that the equation is 4th order.  Therefore, if the
         * root finder returns other than 4 roots, error.
         */
        if ((root_count = rt_poly_roots(C, 4, val)) != 4) 
            break;   /* MISSED */
        

        /* Only real roots indicate an intersection in real space.
         *
         * Look at each root returned; if the imaginary part is zero or
         * sufficiently close, then use the real part as one value of 't'
         * for the intersections
         */
        for (j = 0 ; j < 4; j++) {
        if (NEAR_ZERO(val[j].im, 0.0001)) {
            struct hit *hitp;
            double normalized_dist;
            double distance;
            double3 hit_pt;
            double3 to_hit;
            double angle;

            normalized_dist = val[j].re - cor_proj;
            distance = normalized_dist * bp->bend_radius;

            /* check if this hit is within bend angle */
            hit_pt =  r_pt + distance * r_dir;
            to_hit = hit_pt - vload3(0, bp->bend_V);
            angle = atan2(dot(to_hit, vload3(0, bp->bend_rb)), dot(to_hit, vload3(0, bp->bend_ra)));
            if (angle < 0.0) 
                angle += M_2PI;
            if (angle <= bp->bend_angle) 
            {
            hitp = &hits[*hit_count];
            hitp->hit_dist = distance;
            hitp->hit_vpriv = pprime + normalized_dist * dprime;
            hitp->hit_surfno = seg_no * 10 + PIPE_BEND_OUTER_BODY;

            if ((*hit_count)++ >= RT_PIPE_MAXHITS)
                return;
            }
        }
        }

        if ((bp->bend_alpha_i <= 0.0) || (parallel && dist > bp->bend_ir))
            break;   
        

        /* Now do inner torus */
        A[2] = X2_Y2[2] + cor_pprime.z * cor_pprime.z + 1.0 - bp->bend_alpha_i * bp->bend_alpha_i;

        /* Inline expansion of (void) bn_poly_mul(&Asqr, &A, &A) */
        /* Both polys have degree two */
        Asqr[0] = A[0] * A[0];
        Asqr[1] = A[0] * A[1] + A[1] * A[0];
        Asqr[2] = A[0] * A[2] + A[1] * A[1] + A[2] * A[0];
        Asqr[3] = A[1] * A[2] + A[2] * A[1];
        Asqr[4] = A[2] * A[2];

        /* Inline expansion of bn_poly_scale(&X2_Y2, 4.0) and
         * bn_poly_sub(&C, &Asqr, &X2_Y2).
         */
        C[0] = Asqr[0];
        C[1] = Asqr[1];
        C[2] = Asqr[2] - X2_Y2[0] * 4.0;
        C[3] = Asqr[3] - X2_Y2[1] * 4.0;
        C[4] = Asqr[4] - X2_Y2[2] * 4.0;

        /* It is known that the equation is 4th order.  Therefore,
         * if the root finder returns other than 4 roots, error.
         */
        if ((root_count = rt_poly_roots(C, 4, val)) != 4) 
            break;   /* MISSED */
        

        /* Only real roots indicate an intersection in real space.
         *
         * Look at each root returned; if the imaginary part is zero or
         * sufficiently close, then use the real part as one value of 't'
         * for the intersections
         */
        for (j = 0, root_count = 0; j < 4; j++) {
        if (NEAR_ZERO(val[j].im, 0.0001)) {
            struct hit *hitp;
            double normalized_dist;
            double distance;
            double3 hit_pt;
            double3 to_hit;
            double angle;

            normalized_dist = val[j].re - cor_proj;
            distance = normalized_dist * bp->bend_radius;

            /* check if this hit is within bend angle */
            hit_pt = r_pt + distance * r_dir;
            to_hit = hit_pt - vload3(0, bp->bend_V);
            angle = atan2(dot(to_hit, vload3(0, bp->bend_rb)), dot(to_hit, vload3(0, bp->bend_ra)));
            if (angle < 0.0) 
                angle += M_2PI;
            if (angle <= bp->bend_angle) 
            {
            hitp = &hits[*hit_count];
            hitp->hit_dist = distance;
            hitp->hit_vpriv = pprime + normalized_dist * dprime;
            hitp->hit_surfno = seg_no * 10 + PIPE_BEND_INNER_BODY;

            if ((*hit_count)++ >= RT_PIPE_MAXHITS)
                return;      
            }
        }
        }

        ctr++;
    }


    /* check for surfaces created by discontinuous changes in radii */
    prev = bp->l--;
    if (!prev->pipe_is_bend) {
        struct lin_pipe *lin = (struct lin_pipe *)prev;
        or2_sq = lin->pipe_rotop_sq;
        ir2_sq = lin->pipe_ritop_sq;
        if (!NEAR_EQUAL(or_sq, or2_sq, RT_LEN_TOL) || !NEAR_EQUAL(ir_sq, ir2_sq, RT_LEN_TOL))
           discont_radius_shot(r_pt, r_dir, vload3(0, bp->bend_start), vload3(0, bp->bend_startNorm), or_sq, ir_sq, or2_sq, ir2_sq, hits, hit_count, seg_no);       
    }
    

    next = bp->l++;
    if (next->pipe_is_bend) {
        struct bend_pipe *bend = (struct bend_pipe *)next;
        or2_sq = bend->bend_or * bend->bend_or;
        ir2_sq = bend->bend_ir * bend->bend_ir;
        if (!NEAR_EQUAL(or_sq, or2_sq, RT_LEN_TOL) || !NEAR_EQUAL(ir_sq, ir2_sq, RT_LEN_TOL)) 
           discont_radius_shot(r_pt, r_dir, vload3(0, bp->bend_end), vload3(0, bp->bend_endNorm), or_sq, ir_sq, or2_sq, ir2_sq, hits, hit_count, seg_no);
        
    } else {
        struct lin_pipe *lin = (struct lin_pipe *)next;
        or2_sq = lin->pipe_robase_sq;
        ir2_sq = lin->pipe_ribase_sq;
        if (!NEAR_EQUAL(or_sq, or2_sq, RT_LEN_TOL) || !NEAR_EQUAL(ir_sq, ir2_sq, RT_LEN_TOL)) 
            discont_radius_shot(r_pt, r_dir, vload3(0, bp->bend_end), vload3(0, bp->bend_endNorm), or_sq, ir_sq, or2_sq, ir2_sq, hits, hit_count, seg_no);
        
    }

    return;

}


 void linear_pipe_shot(double3 r_pt, double3 r_dir, struct lin_pipe *lp, struct hit *hits, int *hit_count, int seg_no)
{
    struct hit *hitp;
    double3 work_pt;
    double3 ray_start;
    double3 ray_dir;
    double t_tmp;
    double a, b, c;
    double descrim;

    /* transform ray start point */
    work_pt = r_pt - vload3(0, lp->pipe_V);
    ray_start = MAT4X3VEC(lp->pipe_SoR, work_pt);

    /* rotate ray direction */
    ray_dir = MAT4X3VEC(lp->pipe_SoR, r_dir);

    /* Intersect with outer sides */
    a = ray_dir.x * ray_dir.x + ray_dir.y * ray_dir.y - ray_dir.z * ray_dir.z * lp->pipe_rodiff_sq;
    b = 2.0 * (ray_start.x * ray_dir.x + ray_start.y * ray_dir.y - ray_start.z * ray_dir.z * lp->pipe_rodiff_sq 
        - ray_dir.z * lp->pipe_robase * lp->pipe_rodiff);
    c = ray_start.x * ray_start.x + ray_start.y * ray_start.y - lp->pipe_robase * lp->pipe_robase
    - ray_start.z * ray_start.z * lp->pipe_rodiff_sq - 2.0 * ray_start.z * lp->pipe_robase * lp->pipe_rodiff;

    descrim = b * b - 4.0 * a * c;

    if (descrim > 0.0) {
    double sqrt_descrim;
    double3 hit_pt;

    sqrt_descrim = sqrt(descrim);

    t_tmp = (-b - sqrt_descrim) / (2.0 * a);
    hit_pt = ray_start + t_tmp * ray_dir;
    if (hit_pt.z >= 0.0 && hit_pt.z <= 1.0) 
    {
        hitp = &hits[*hit_count];
        hitp->hit_dist = t_tmp;
        hitp->hit_surfno = seg_no * 10 + PIPE_LINEAR_OUTER_BODY;
        hitp->hit_vpriv = vload3(0, hit_pt);
        hitp->hit_vpriv.z = (-lp->pipe_robase - hit_pt.z * lp->pipe_rodiff) * lp->pipe_rodiff;

        if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
            return;
        
    }

    t_tmp = (-b + sqrt_descrim) / (2.0 * a);
    hit_pt = ray_start + t_tmp * ray_dir;
    if (hit_pt.z >= 0.0 && hit_pt.z <= 1.0) 
    {
        hitp = &hits[*hit_count];
        hitp->hit_dist = t_tmp;
        hitp->hit_surfno = seg_no * 10 + PIPE_LINEAR_OUTER_BODY;
        hitp->hit_vpriv = vload3(0, hit_pt);
        hitp->hit_vpriv.z = (-lp->pipe_robase - hit_pt.z * lp->pipe_rodiff) * lp->pipe_rodiff;

        if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
            return;
        
    }
    }

    if (lp->pipe_ribase > 0.0 || lp->pipe_ritop > 0.0) {
    /* Intersect with inner sides */

    a = ray_dir.x * ray_dir.x + ray_dir.y * ray_dir.y - ray_dir.z * ray_dir.z * lp->pipe_ridiff_sq;
    b = 2.0 * (ray_start.x * ray_dir.x + ray_start.y * ray_dir.y - ray_start.z * ray_dir.z * lp->pipe_ridiff_sq
           - ray_dir.z * lp->pipe_ribase * lp->pipe_ridiff);
    c = ray_start.x * ray_start.x + ray_start.y * ray_start.y - lp->pipe_ribase * lp->pipe_ribase
        - ray_start.z * ray_start.z * lp->pipe_ridiff_sq - 2.0 * ray_start.z * lp->pipe_ribase * lp->pipe_ridiff;

    descrim = b * b - 4.0 * a * c;

    if (descrim > 0.0) {
        double sqrt_descrim;
        double3 hit_pt;

        sqrt_descrim = sqrt(descrim);

        t_tmp = (-b - sqrt_descrim) / (2.0 * a);
        hit_pt = ray_start + t_tmp * ray_dir;
        if (hit_pt.z >= 0.0 && hit_pt.z <= 1.0) {
        hitp = &hits[*hit_count];
        hitp->hit_dist = t_tmp;
        hitp->hit_surfno = seg_no * 10 + PIPE_LINEAR_INNER_BODY;
        hitp->hit_vpriv = vload3(0, hit_pt);
        hitp->hit_vpriv.z = (-lp->pipe_ribase - hit_pt.z * lp->pipe_ridiff) * lp->pipe_ridiff;

        if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
            return;
        
        }

        t_tmp = (-b + sqrt_descrim) / (2.0 * a);
        hit_pt = ray_start + t_tmp * ray_dir;
        if (hit_pt.z >= 0.0 && hit_pt.z <= 1.0) {
        hitp = &hits[*hit_count];
        hitp->hit_dist = t_tmp;
        hitp->hit_surfno = seg_no * 10 + PIPE_LINEAR_INNER_BODY;
        hitp->hit_vpriv = vload3(0, hit_pt);
        hitp->hit_vpriv.z = (-lp->pipe_ribase - hit_pt.z * lp->pipe_ridiff) * lp->pipe_ridiff;

        if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
            return;
        
        }
    }
    }

}







 void
pipe_start_shot(const double3 r_pt, const double3 r_dir, struct id_pipe *id_p, struct hit *hits, int *hit_count, int seg_no)
{
    double3 hit_pt;
    double t_tmp;
    double radius_sq;

    if (!id_p->pipe_is_bend) 
    {
        struct lin_pipe *lin = (struct lin_pipe *)(&id_p->linear);
        double dist_to_plane;
        double norm_dist;
        double slant_factor;

        dist_to_plane = dot(vload3(0, lin->pipe_H), vload3(0, lin->pipe_V));
        norm_dist = dist_to_plane - dot(vload3(0, lin->pipe_H), r_pt);
        slant_factor = dot(vload3(0, lin->pipe_H), r_dir);
        if (!ZERO(slant_factor)) 
        {
            double3 to_center;

            t_tmp = norm_dist / slant_factor;
            hit_pt = r_pt + t_tmp * r_dir;
            to_center = vload3(0, lin->pipe_V) - hit_pt;
            radius_sq = dot(to_center, to_center);
            if (radius_sq <= lin->pipe_robase_sq && radius_sq >= lin->pipe_ribase_sq) 
            {
            hits[*hit_count].hit_dist = t_tmp;
            hits[*hit_count].hit_surfno = seg_no * 10 + PIPE_LINEAR_BASE;

            if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
                return;
            }
        }
    } else if (id_p->pipe_is_bend) 
    {
        struct bend_pipe *bend = (struct bend_pipe *)(&id_p->curved);
        double dist_to_plane;
        double norm_dist;
        double slant_factor;

        dist_to_plane = dot(vload3(0, bend->bend_rb), vload3(0, bend->bend_start));
        norm_dist = dist_to_plane - dot(vload3(0, bend->bend_rb), r_pt);
        slant_factor = dot(vload3(0, bend->bend_rb), r_dir);

        if (!ZERO(slant_factor)) 
        {
            double3 to_center;

            t_tmp = norm_dist / slant_factor;
            hit_pt = r_pt + t_tmp * r_dir;
            to_center = vload3(0, bend->bend_start) - hit_pt;
            radius_sq = dot(to_center, to_center);
            if (radius_sq <= bend->bend_or * bend->bend_or && radius_sq >= bend->bend_ir * bend->bend_ir) {
            
            hits[*hit_count].hit_dist = t_tmp;
            hits[*hit_count].hit_surfno = seg_no * 10 + PIPE_BEND_BASE;

            if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
                return;
            
            }
        }
    }
}






 void
pipe_end_shot(const double3 r_pt, const double3 r_dir, struct id_pipe *id_p, struct hit *hits, int *hit_count, int seg_no)
{
    double3 hit_pt;
    double t_tmp;
    double radius_sq;
    struct hit *hitp;

    if (!id_p->pipe_is_bend) {
    struct lin_pipe *lin = (struct lin_pipe *)(&id_p->linear);
    double3 top;
    double dist_to_plane;
    double norm_dist;
    double slant_factor;

    top = vload3(0, lin->pipe_V) + lin->pipe_len * vload3(0, lin->pipe_H);
    dist_to_plane = dot(vload3(0, lin->pipe_H), top);
    norm_dist = dist_to_plane - dot(vload3(0, lin->pipe_H), r_pt);
    slant_factor = dot(vload3(0, lin->pipe_H), r_dir);
    if (!ZERO(slant_factor)) {
        double3 to_center;

        t_tmp = norm_dist / slant_factor;
        hit_pt = r_pt + t_tmp * r_dir;
        to_center = top - hit_pt;
        radius_sq = dot(to_center, to_center);
        
        if (radius_sq <= lin->pipe_rotop_sq && radius_sq >= lin->pipe_ritop_sq) {
        hits[*hit_count].hit_dist = t_tmp;
        hits[*hit_count].hit_surfno = seg_no * 10 + PIPE_LINEAR_TOP;

        if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
            return;
        }
    }
    } else if (id_p->pipe_is_bend) {
    struct bend_pipe *bend = (struct bend_pipe *)(&id_p->curved);
    double3 to_end;
    double3 plane_norm;
    double dist_to_plane;
    double norm_dist;
    double slant_factor;

    to_end = vload3(0, bend->bend_end) - vload3(0, bend->bend_V);
    plane_norm = cross(to_end, vload3(0, bend->bend_N));
    plane_norm = normalize(plane_norm);

    dist_to_plane = dot(plane_norm, vload3(0,bend->bend_end));
    norm_dist = dist_to_plane - dot(plane_norm, r_pt);
    slant_factor = dot(plane_norm, r_dir);

    if (!ZERO(slant_factor)) {
        double3 to_center;

        t_tmp = norm_dist / slant_factor;
        hit_pt = r_pt + t_tmp * r_dir;
        to_center = vload3(0, bend->bend_end) - hit_pt;
        radius_sq = dot(to_center, to_center);
        if (radius_sq <= bend->bend_or * bend->bend_or && radius_sq >= bend->bend_ir * bend->bend_ir) {
        hits[*hit_count].hit_dist = t_tmp;
        hits[*hit_count].hit_surfno = seg_no * 10 + PIPE_BEND_TOP;

        if ((*hit_count)++ >= RT_PIPE_MAXHITS) 
            return;
        }
    }
    }
}



void rt_pipe_elim_dups(struct hit *hit, int *nh, const double3 r_pt, const double3 r_dir)
{
    struct hit *hitp;
    struct hit *next_hit;
    int hitNo = 0;

    /* delete duplicate hits */
    while (hitNo < ((*nh) - 1)) 
    {
        hitp = &hit[hitNo];
        next_hit = &hit[hitNo + 1];

        if (NEAR_EQUAL(hitp->hit_dist, next_hit->hit_dist, RT_PCOEF_TOL) && hitp->hit_surfno == next_hit->hit_surfno) 
        {
            for (int i = hitNo ; i < (*nh-2) ; i++) 
                hit[i] = hit[i + 2];
            
            (*nh)=(*nh)-2;
        } 
        else 
            hitNo++;
    }

    if ((*nh) == 1) {
    (*nh) = 0;
    return;
    }

    if ((*nh) == 0 || (*nh) == 2) 
    return;

    hitNo = 0;
    while (hitNo < ((*nh) - 1)) 
    {
        int hitNoPlus = hitNo + 1;
        /* struct hit *first = &hit[hitNo]; */
        struct hit *second = &hit[hitNoPlus];

        /* keep first entrance hit, eliminate all successive entrance hits */
        while (hitNoPlus < (*nh) && dot(second->hit_normal, r_dir) < 0.0) 
        {
            for (int j = hitNoPlus ; j < ((*nh) - 1) ; j++)
                hit[j] = hit[j + 1];
            
            (*nh)--;
            second = &hit[hitNoPlus];
        }

        /* second is now an exit hit at hit[hitNoPlus] */

        /* move to next hit */
        hitNoPlus++;
        if (hitNoPlus >= (*nh))
            break;
        
        /* set second to the next hit */
        second = &hit[hitNoPlus];

        /* eliminate all exit hits (except the last one) till we find another entrance hit */
        while (hitNoPlus < (*nh) && dot(second->hit_normal, r_dir) > 0.0) {
            int j;
            for (j = hitNoPlus - 1 ; j < ((*nh) - 1) ; j++) {
            hit[j] = hit[j + 1];
            }
            (*nh)--;
            second = &hit[hitNoPlus];
        }
        hitNo = hitNoPlus;
    }
}



/**
 * Given ONE ray distance, return the normal and entry/exit point.
 */
pipe_norm(struct hit *hitp, const double3 r_pt, const double3 r_dir, global const struct pipe_specific *pipe)
{
    struct id_pipe *pipe_id;
    pipe_id = pipe->pipe_id;

    struct lin_pipe *pipe_lin;
    struct bend_pipe *pipe_bend;
    double w;
    double3 work;
    double3 work1;

    int segno = hitp->hit_surfno / 10;

    for (int i = 1; i < segno; i++) 
        pipe_id = pipe_id->l++;
    

    pipe_lin = (struct lin_pipe *)pipe_id;
    pipe_bend = (struct bend_pipe *)pipe_id;

    hitp->hit_point = r_pt + hitp->hit_dist * r_dir;

    switch (hitp->hit_surfno % 10) {
    case PIPE_LINEAR_TOP:
        hitp->hit_normal = vload3(0, pipe_lin->pipe_H);
        break;

    case PIPE_LINEAR_BASE:
        hitp->hit_normal = -vload3(0, pipe_lin->pipe_H);
        break;

    case PIPE_LINEAR_OUTER_BODY:
        hitp->hit_normal = MAT4X3VEC(pipe_lin->pipe_invRoS, hitp->hit_vpriv);
        hitp->hit_normal = normalize(hitp->hit_normal);
        break;

    case PIPE_LINEAR_INNER_BODY:
        hitp->hit_normal = MAT4X3VEC(pipe_lin->pipe_invRoS, hitp->hit_vpriv);
        hitp->hit_normal = normalize(hitp->hit_normal);
        hitp->hit_normal = -vload3(0, hitp->hit_normal);
        break;
    case PIPE_BEND_OUTER_BODY:
        w = hitp->hit_vpriv.x * hitp->hit_vpriv.x + hitp->hit_vpriv.y * hitp->hit_vpriv.y + hitp->hit_vpriv.z * hitp->hit_vpriv.z +
        1.0 - pipe_bend->bend_alpha_o * pipe_bend->bend_alpha_o;
        work = (double3)((w - 2.0) * hitp->hit_vpriv.x, (w - 2.0) * hitp->hit_vpriv.y, w * hitp->hit_vpriv.z);
        work = normalize(work);
        hitp->hit_normal = MAT3X3VEC(pipe_bend->bend_invR, work);
        break;
    case PIPE_BEND_INNER_BODY:
        w = hitp->hit_vpriv.x * hitp->hit_vpriv.x + hitp->hit_vpriv.y * hitp->hit_vpriv.y + hitp->hit_vpriv.z * hitp->hit_vpriv.z +
        1.0 - pipe_bend->bend_alpha_i * pipe_bend->bend_alpha_i;
        work = (double3) ((w - 2.0) * hitp->hit_vpriv.x, (w - 2.0) * hitp->hit_vpriv.y, w * hitp->hit_vpriv.z);
        work = normalize(work);
        work1 = MAT3X3VEC(pipe_bend->bend_invR, work);
        hitp->hit_normal = - vload3(0, work1);
        break;
    case PIPE_BEND_BASE:
        hitp->hit_normal = - vload3(0, pipe_bend->bend_rb);
        break;
    case PIPE_BEND_TOP:
        work = pipe_bend->bend_end - pipe_bend->bend_V;
        hitp->hit_normal = cross(pipe_bend->bend_N, work);
        hitp->hit_normal = normalize(hitp->hit_normal);
        break;
    case PIPE_RADIUS_CHANGE:
        break; /* already have normal */
    default:
        break;
    }
}






/**
 * Intersect a ray with a pipe.  If an intersection occurs, a struct
 * seg will be acquired and filled in.
 *
 * Returns -
 *  0 MISS
 * >0 HIT
 */
int pipe_shot(RESULT_TYPE *res, const double3 r_pt, const double3 r_dir, const uint idx, global const struct pipe_specific *pipe)
{
    struct id_pipe *pipe_id;
    pipe_id = pipe->pipe_id;
    struct hit hits[RT_PIPE_MAXHITS];
    int total_hits = 0;
    int seg_no = pipe->seg_no;
    double3 invdir;
    int i;

    /* compute the inverse of the direction cosines */
    if (!ZERO(r_dir.x)) 
        invdir.x = 1.0/r_dir.x;
    else 
        invdir.x = INFINITY;
    if (!ZERO(r_dir.y)) 
        invdir.y = 1.0/r_dir.y;
    else
        invdir.y = INFINITY;
    if (!ZERO(r_dir.z)) 
        invdir.z = 1.0/r_dir.z;
    else 
        invdir.z = INFINITY;

    pipe_start_shot(r_pt, r_dir, pipe_id, hits, &total_hits, 1);
    pipe_end_shot(r_pt, r_dir, pipe_id+seg_no, hits, &total_hits, seg_no);

    for (i=0;i <=seg_no;i++, pipe_id++) 
    {
    	if (!pipe_id->pipe_is_bend) {
    	    struct lin_pipe *lin = (struct lin_pipe *)pipe_id;
    	    if (!rt_in_rpp(r_pt, invdir, lin->pipe_min, lin->pipe_max)) 
    		  continue;
    	    
    	    linear_pipe_shot(r_pt, r_dir, lin, hits, &total_hits, i+1);
    	} else {
    	    struct bend_pipe *bend = (struct bend_pipe *)pipe_id;
    	    if (!rt_in_sph(r_pt, r_dir, bend->bend_bound_center, bend->bend_bound_radius_sq))
    		  continue;
    	    
    	    bend_pipe_shot(r_pt, r_dir, bend, hits, &total_hits, i+1);
    	}
    }

    if (!total_hits)
	   return 0;

    /* calculate hit points and normals */
    for (i = 0 ; i < total_hits ; i++) 
	   pipe_norm(&hits[i], r_pt, r_dir, pipe);

    /* sort the hits */
    primitive_hitsort(hits, total_hits);

    /* eliminate duplicate hits */
    rt_pipe_elim_dups(hits, &total_hits, r_pt, r_dir);

    /* Bad number of hits */
    if (total_hits % 2) 
	   return 0;

    for (i = 0 ; i < total_hits ; i += 2) 
	   do_segp(res, idx, &hits[i], &hits[i+1]);

    if (total_hits) 
	   return 1;    /* HIT */
    else 
	   return 0;    /* MISS */
    
}




/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */