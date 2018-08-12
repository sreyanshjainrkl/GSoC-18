#include "common.cl"

#define XMIN 0
#define XMAX 1
#define YMIN 2
#define YMAX 3
#define ZMIN 4
#define ZMAX 5
#define ZMID 6
#define ZTOP 7
#define DSP_CUT_DIR_ADAPT 'a'
#define DSP_CUT_DIR_llUR 'l'
#define DSP_CUT_DIR_ULlr 'L'


static const double3 dsp_pl[7] = {
    {-1.0, 0.0, 0.0},
    { 1.0, 0.0, 0.0},

    {0.0, -1.0, 0.0},
    {0.0,  1.0, 0.0},

    {0.0, 0.0, -1.0},
    {0.0, 0.0,  1.0},
    {0.0, 0.0,  1.0},
};

struct dsp_bb {
    int dsp_min[3];
    int dsp_max[3];
    int dspb_subcell_size;/* XXX This is not yet computed */
    int dspb_ch_dim[2];	/* dimensions of children[] */
    struct dsp_bb *dspb_children[16];
};

struct dsp_specific {
    int xsiz;
    int ysiz;
    struct dsp_bb *bb_array;
    int *dsp_buf;
    double dsp_mtos[16];  
    double dsp_stom[16];
    char dsp_cuttype;
    int dsp_smooth;
    int	dsp_xcnt;
    int	dsp_ycnt;
};

#define DSP(_p, _x, _y) (((_p) && (_p)->dsp_buf) ?	((int *)((_p)->dsp_buf))[    \
	    (_y) * ((struct dsp_specific *)_p)->dsp_xcnt + (_x)] : 0)


static inline int
dsp_segp(RESULT_TYPE *res, const uint idx, struct hit *in_hit, struct hit *out_hit, struct dsp_specific *dsp, double3 r_dir)
{
	double3 r_dir1, dir, v;

    // if both points are on the "floor" of the DSP, then we don't have a hit segment	
	if (NEAR_ZERO(in_hit->hit_point.z, RT_PCOEF_TOL) && NEAR_ZERO(out_hit->hit_point.z, RT_PCOEF_TOL))
		return 0;
		    
	// throw away any zero length segments, mostly to avoid seeing inside-out segments
	if (ZERO(out_hit->hit_dist - in_hit->hit_dist))
		return 0;

    if (in_hit->hit_surfno < ZMAX) {
		in_hit->hit_vpriv.x = 0.0;
		in_hit->hit_vpriv.y = 0.0;
    }

    if (out_hit->hit_surfno < ZMAX) {
		out_hit->hit_vpriv.x = 0.0;
		out_hit->hit_vpriv.y = 0.0;
    }

    in_hit->hit_vpriv.z = 0.0; /* flag as in-hit */
    out_hit->hit_vpriv.z = 1.0; /* flag as out-hit */

    r_dir1 = vload3(0, MAT4X3VEC(dsp->dsp_mtos, r_dir));
    r_dir1 = normalize(r_dir1);

	dir = r_dir1 * in_hit->hit_dist;
	v = MAT4X3VEC(dsp->dsp_stom, dir);
	in_hit->hit_dist = length(v);
	if (dot(v, r_dir) < 0.0) in_hit->hit_dist *= -1.0;

	dir = r_dir1 * out_hit->hit_dist;
	v = MAT4X3VEC(dsp->dsp_stom, dir);
	out_hit->hit_dist = length(v);
	if (dot(v, r_dir) < 0.0) out_hit->hit_dist *= -1.0;
	    
    do_segp(res, idx, &in_hit, &out_hit);		
    return 2;
}



int dsp_in_rpp(double3 pt, double3 invdir, double3 min, double3 max, double *rmin, double *rmax, int *dmin, int *dmax)
{

    *rmin = -MAX_FASTF;
    *rmax =  MAX_FASTF;
    *dmin = -1;
    *dmax = -1;

    /* Start with infinite ray, and trim it down */

    double x0 = (min.x - pt.x) * invdir.x;
    double y0 = (min.y - pt.y) * invdir.y;
    double z0 = (min.z - pt.z) * invdir.z;
    double x1 = (max.x - pt.x) * invdir.x;
    double y1 = (max.y - pt.y) * invdir.y;
    double z1 = (max.z - pt.z) * invdir.z;


    /* X axis */
    if (invdir.x < 0.0) {
        if (*rmax > x0) {
            *rmax = x0;
            *dmax = XMIN;
        }
        if (*rmin < x1) {
            *rmin = x1;
            *dmin = XMAX;
        }
    }  else if (invdir.x > 0.0) {
        if (*rmax > x1) {
            *rmax = x1;
            *dmax = XMAX;
        }
        if (*rmin < x0) {
            *rmin = x0;
            *dmin = XMIN;
        }
    } else {
        if ((min.x > pt.x) || (max.x < pt.x))
            return 0;   /* MISS */
    }

    /* Y axis */
    if (invdir.y < 0.0) {
    if (*rmax > y0) {
        /* towards smaller */
        *rmax = y0;
        *dmax = YMIN;
    }
    if (*rmin < y1) {
        *rmin = y1;
        *dmin = YMAX;
    }
    }  else if (invdir.y > 0.0) {
    /* towards larger */
    if (*rmax > y1) {
        *rmax = y1;
        *dmax = YMAX;
    }
    if (*rmin < y0) {
        *rmin = y0;
        *dmin = YMIN;
    }
    } else {
        if ((min.y > pt.y) || (max.y < pt.y))
            return 0;   /* MISS */
    }

    /* Z axis */
    
    if (invdir.z < 0.0) {
    /* towards smaller */
    if (*rmax > z0) {
        *rmax = z0;
        *dmax = ZMIN;
    }
    if (*rmin < z1) {
        *rmin = z1;
        *dmin = ZMAX;
    }
    }  else if (invdir.z > 0.0) {
    /* towards larger */
    if (*rmax > z1) {
        *rmax = z1;
        *dmax = ZMAX;
    }
    if (*rmin < z0) {
        *rmin = z0;
        *dmin = ZMIN;
    }
    } else {
        if ((min.z > pt.z) || (max.z < pt.z))
            return 0;   /* MISS */
    }

    /* If equal, RPP is actually a plane */
    if (*rmin > *rmax)
       return 0;   /* MISS */
    else
       return 1;       /* HIT */
}




void permute_cell(double3 *A, double3 *B, double3 *C, double3 *D, struct dsp_specific *dsp, struct dsp_bb *dsp_bb)
{
    
    switch (dsp->dsp_cuttype) 
    {
	case DSP_CUT_DIR_llUR:
	    break;

	case DSP_CUT_DIR_ADAPT: 
	 {
	    int lo[2], hi[2];
	    double h1, h2, h3, h4;
	    double cAD, cBC;  /* curvature in direction AD, and BC */
	    double3 *tmp;

	    lo[X] = fmax(dsp_bb->dsp_min[X] - 1, 0);
	    lo[Y] = fmax(dsp_bb->dsp_min[Y] - 1, 0);
	    hi[X] = fmin(dsp_bb->dsp_max[X] + 1, dsp->xsiz);
	    hi[Y] = fmin(dsp_bb->dsp_max[Y] + 1, dsp->ysiz);

	    /* compute curvature along the A->D direction */
	    h1 = DSP(dsp, lo[X], lo[Y]);
	    h2 = A.z;
	    h3 = D.z;
	    h4 = DSP(dsp, hi[X], hi[Y]);

	    cAD = fabs(h3 + h1 - 2.0*h2) + fabs(h4 + h2 - 2.0*h3);

	    /* compute curvature along the B->C direction */
	    h1 = DSP(dsp, hi[X], lo[Y]);
	    h2 = B.z;
	    h3 = C.z;
	    h4 = DSP(dsp, lo[X], hi[Y]);

	    cBC = fabs(h3 + h1 - 2.0*h2) + fabs(h4 + h2 - 2.0*h3);

	    if (cAD < cBC) 
	    	break;

	    /* prefer the B-C cut */
	    tmp = A;
	    A = B;
	    B = D;
	    D = C;
	    C = tmp;
	    break;
	 }
	case DSP_CUT_DIR_ULlr:
	 {
	 	double3 *tmp;

	 	tmp = A;
	 	A = B;
	 	B = tmp;
	 	tmp = C;
	 	C = D;
	 	D = tmp;
	    break;
     }
    }
}



int isect_ray_triangle(double3 r_pt1, double3 r_dir1, double3 A, double3 B, double3 C, struct hit *hitp)
{
    double3 P;			/* plane intercept point */
    double3 AB, AC, AP;
    double3 N;			/* Normal for plane of triangle */
    double NdotDir;
    double alpha, beta;	/* barycentric distances */
    double hitdist;		/* distance to ray/triangle intercept */

    AB = B - vload3(0, A);
    AC = C - vload3(0, A);

    /* Compute the plane equation of the triangle */
    N = cross(AB, AC);
    N = normalize(N);

    /* intersect ray with plane */
    NdotDir = dot(N, r_dir1);
    /* Ray perpendicular to plane of triangle */
    if (fabs(NdotDir)<(1e-6)) 
		return -1;

    /* dist to plane icept */
    hitdist = (dot(N, A)-dot(N, r_pt1)) / NdotDir;

    P = r_pt1 + hitdist * r_dir1;
    AP = P - vload3(A);

    if (ZERO(AB.x)) 
		beta = AB.y * AP.y;
    else 
		beta = AB.x * AP.x;
    
    if (ZERO(AC.x))
		alpha = AC.y * AP.y;
    else 
		alpha = AC.x * AP.x;

    if (alpha < -RT_PCOEF_TOL || beta < -RT_PCOEF_TOL || (alpha+beta) > (1.0 + RT_PCOEF_TOL))
		return 0;

    hitp->hit_dist = hitdist;
    hitp->hit_normal = vload3(N);
    hitp->hit_point = vload3(P);
    return 1;
}


int check_bbpt_hit_elev(int i, double3 A, double3 B, double3 C, double3 D, double3 P)
{
    double slope = 0.0;
    double delta = 0.0;
    double origin = 0.0;

    switch (i) {
	
	case XMIN:
	    slope = C.z - A.z;
	    delta = P.y - A.y;
	    origin = A.z;
	    break;

	case XMAX:
	    slope = D.z - B.z;
	    delta = P.y - B.y;
	    origin = B.z;
	    break;
	
	case YMIN:
	    slope = B.z - A.z;
	    delta = P.x - A.x;
	    origin = A.z;
	    break;
	
	case YMAX:
	    slope = D.z - C.z;
	    delta = P.x - C.x;
	    origin = C.z;
	    break;
	
	case ZMIN:
	    return 1;
	    break;
	
	case ZMAX:
	    return 0;
	    break;
	
	default:
	    break;
    }

    if ((origin + slope * delta) < P.z) return 0;

    return 1;
}



int dsp_shot(RESULT_TYPE *res, const double3 r_pt, const double3 r_dir, const uint idx, global const struct dsp_specific *dsp)
{
    double3 r_pt1, r_dir1; 
    double3 invdir;
    struct hit hits[4];
    int dmin, dmax;
    double r_min, r_max;
    int a_onehit=0;

    double3 bbmin, bbmax;
    double3 minpt, maxpt;
    int min_z;
    
    int nhits=0;
    int npts=0;



    /* map ray into the coordinate system of the dsp  */
    r_pt1 = MAT4X3PNT(vload16(0, dsp->dsp_mtos), r_pt);
    r_dir1 = MAT4X3VEC(dsp->dsp_mtos, r_dir);
    r_dir1 = normalize(r_dir1);


    /* compute the inverse of the direction cosines */
    if (!ZERO(r_dir1.x)) 
		invdir.x = 1.0/r_dir1.x;
    else 
		invdir.x = MAX_FASTF;

    if (!ZERO(r_dir1.y)) 
		invdir.y = 1.0/r_dir1.y;
    else 
		invdir.y = MAX_FASTF;
    
    if (!ZERO(r_dir1.z)) 
		invdir.z = 1.0/r_dir1.z;
    else 
		invdir.z = MAX_FASTF;
    
	int ctr=0;
	struct dsp_bb **nodes;
	nodes = &dsp->bb_array;
	ctr++;

	while(ctr>0)
	{
		ctr--;
		struct dsp_bb *node;
		node = *(nodes+ctr); 

	    /* check to see if we miss the RPP for this area entirely */
	    bbmax = vload3(0, node->dsp_max);
	    bbmin = double3(node->dsp_min[X], node->dsp_min[Y], 0.0);

	    if (! dsp_in_rpp(r_pt1, invdir, bbmin, bbmax, &r_min, &r_max, &dmin, &dmax))
			break;    /* missed it all, just return */
	 

	    /* At this point we know that we've hit the overall bounding box  */

	    minpt = r_pt1 + r_min * r_dir1;
	    maxpt = r_pt1 + r_max * r_dir1;

	    /* if both hits are UNDER the top of the "foundation" pillar, we
	     * can just add a segment for that range and return
	     */
	    min_z = node->dsp_min[Z];

	    if (minpt.z < min_z && maxpt.z < min_z) 
	    {
			hits[0].hit_vpriv = (double3) (0.0);
			hits[1].hit_vpriv = (double3) (0.0);

			hits[0].hit_dist = r_min;
			hits[0].hit_point = vload3(0, minpt);
			hits[0].hit_normal = (double3)(dsp_pl[dmin]);
			hits[0].hit_surfno = dmin;		

			hits[1].hit_dist = r_max;
			hits[1].hit_point = vload3(0, maxpt);
			hits[1].hit_normal = (double3)(dsp_pl[dmax]);
			hits[1].hit_surfno = dmax;

			npts+= dsp_segp(res, idx, &hits[0], &hits[1], dsp, r_dir);
			nhits++;
			if (( r_min > 0.0 || r_max> 0.0 )&&(nhits > a_onehit))
				return npts;
		
	    }

	    /* We've hit something where we might be going through the
	    * boundary.  We've got to intersect the children
	    */
		else if (node->dspb_ch_dim[0]) 
			{

			    double tDX;		/* dist along ray to span 1 cell in X dir */
			    double tDY;		/* dist along ray to span 1 cell in Y dir */
			    double tX, tY;	/* dist from hit pt. to next cell boundary */
			    double node_dist;
			    short cX, cY;	/* coordinates of current cell */
			    short cs;		/* cell X, Y dimension */
			    short stepX, stepY;	/* dist to step in child array for each dir */
			    short stepPX, stepPY;
			    double out_dist;
			    struct dsp_bb **p;
			 
			    /* compute the size of a cell in each direction */
			    cs = node->dspb_subcell_size;

			    /* compute current cell */
			    cX = (minpt.x - bbmin.x) / cs;
			    cY = (minpt.y - bbmin.y) / cs;

			    /* a little bounds checking because a hit on XMAX or YMAX looks
			     * like it should be in the next cell outside the box
			     */
			    if (cX >= node->dspb_ch_dim[X]) cX = node->dspb_ch_dim[X] - 1;
			    if (cY >= node->dspb_ch_dim[Y]) cY = node->dspb_ch_dim[Y] - 1;


			    tX = tY = node_dist = r_min;

			    if (r_dir1.x < 0.0) {
				stepPX = stepX = -1;
				tDX = -cs / r_dir1.x;
				tX += ((bbmin.x + (cX * cs)) - minpt.x) / r_dir1.x;
			    } 

			    else {
				stepPX = stepX = 1;
				tDX = cs / r_dir1.x;

				if (r_dir1.x > 0.0)
				    tX += ((bbmin.x + ((cX+1) * cs)) - minpt.x) / r_dir1.x;
				else
				    tX = MAX_FASTF; /* infinite distance to next X boundary */
			    }

			    if (r_dir1.y < 0) {
				stepY = -1;
				stepPY = -node->dspb_ch_dim[X];
				tDY = -cs / r_dir1.y;
				tY += ((bbmin.y + (cY * cs)) - minpt.y) / r_dir1.y;
			    } else {
				stepY = 1;
				stepPY = node->dspb_ch_dim[X];
				tDY = cs / r_dir1.y;

				if (r_dir1.y > 0.0)
				    tY += ((bbmin.y + ((cY+1) * cs)) - minpt.y) / r_dir1.y;
				else
				    tY = MAX_FASTF;
			    }

			    /* factor in the tolerance to the out-distance */
			    out_dist = r_max - RT_PCOEF_TOL;

			    p = &node->dspb_children[node->dspb_ch_dim[X] * cY + cX];

			    do 
			    {
					*(nodes+ctr) = *p;
					ctr++;

					/* figure out which cell is next */
					if (tX < tY) 
						{
						    cX += stepX;  /* track cell offset for debugging */
						    p += stepPX;

						    node_dist = tX;
						    tX += tDX;
						} 
					else 
						{
						    cY += stepY;  /* track cell offset for debugging */
						    p += stepPY;
						    node_dist = tY;
						    tY += tDY;
						}

			    } while (node_dist < out_dist && cX < node->dspb_ch_dim[X] && cX >= 0 && cY < node->dspb_ch_dim[Y] && cY >= 0);

			}
		    else
			{
			    /* This section is for level 0 intersections only  */

			    bbmin.z = node->dsp_min[Z];
			    if (dsp_in_rpp(r_pt1, invdir, bbmin, bbmax, &r_min, &r_max, &dmin, &dmax)) 
					{
						double3 A, B, C, D, P;
					    int x, y;
					    struct hit *hitp;
					    int hitf = 0;	/* bit flags for valid hits in hits */
					    int cond, i;
					    int hitcount = 0;
					    double dot1, dot2;


					    x = node->dsp_min[X];
					    y = node->dsp_min[Y];
					    A = double3(x, y, DSP(dsp, x, y));

					    x = node->dsp_max[X];
					    B = double3(x, y, DSP(dsp, x, y));

					    y = node->dsp_max[Y];
					    D = double3(x, y, DSP(dsp, x, y));

					    x = node->dsp_min[X];
					    C = double3(x, y, DSP(dsp, x, y));



					    /* first order of business is to discard any "fake" hits on the
					     * bounding box, and fill in any "real" hits in our list
					     */
					    P = r_pt1 + r_min * r_dir1;
					    if (check_bbpt_hit_elev(dmin, A, B, C, D, P)) 
					    {
							hits[0].hit_dist = r_min;
							hits[0].hit_point = vload3(0, P);
							hits[0].hit_normal = vload3(0, dsp_pl[dmin]);
							/* vpriv */
							hits[0].hit_vpriv.x = node->dsp_min[X];
							hits[0].hit_vpriv.y = node->dsp_min[Y];
							/* private */
							hits[0].hit_surfno = dmin;

							hitcount++;
							hitf = 1;
						}

					    /* make sure the point P is below the cell top */
					    P = r_pt1 + r_max * r_dir1;
					    if (check_bbpt_hit_elev(dmax, A, B, C, D, P)) 
					    {
							/* P is at or below the top surface */
							hits[3].hit_dist = r_max;
							hits[3].hit_point = vload3(P);
							hits[3].hit_normal = vload3(dsp_pl[dmax]);
							/* vpriv */
							hits[3].hit_vpriv.x = node->dsp_min[X];
							hits[3].hit_vpriv.y = node->dsp_min[Y];
							/* private */
							hits[3].hit_surfno = dmax;

							hitcount++;
							hitf |= 8;
						}

						permute_cell(&A, &B, &C, &D, dsp, node);

					    if ((cond=isect_ray_triangle(r_pt1, r_dir1, B, D, A, &hits[1])) > 0.0) 
					    {
							hits[1].hit_vpriv.x = node->dsp_min[X];
							hits[1].hit_vpriv.y = node->dsp_min[Y];
							hits[1].hit_surfno = ZTOP; /* indicate we hit the top */

							hitcount++;
							hitf |= 2;
					    } 

					    if ((cond=isect_ray_triangle(r_pt1, r_dir1, C, A, D, &hits[2])) > 0.0) 
					    {
							hits[2].hit_vpriv.x = node->dsp_min[X];
							hits[2].hit_vpriv.y = node->dsp_min[Y];
							hits[2].hit_surfno = ZTOP; /* indicate we hit the top */

							hitcount++;
							hitf |= 4;

							if (hitf & 2) 
							{
							    /* if this hit occurs before the hit on the other triangle
							     * swap the order
							     */
							    if (hits[1].hit_dist > hits[2].hit_dist) {
								struct hit *tmp;
								tmp = hits[1]; /* struct copy */
								hits[1] = hits[2]; /* struct copy */
								hits[2] = tmp; /* struct copy */

							    } 
							}


					    } 

					    /* fill out the segment structures */

					    hitp = 0;
					    for (i = 0; i < 4; i++) 
					    {
						if (hitf & (1<<i)) 
						{
						    if (hitp) 
						    {
							dot2 = dot(r_dir1, hits[i].hit_normal);
							/* if we have two entry points then pick the first one */
							if (dot2 < 0.0) {
							    if (hitp->hit_dist > hits[i].hit_dist)
								     hitp = &hits[i];
								continue;
							}

							/* create seg with hits[i].hit_point as out point */
							npts+ = dsp_segp(res, idx, hitp, &hits[i], dsp, r_dir)
							nhits++;
							if (( r_min > 0.0 || r_max> 0.0 )&&(nhits > a_onehit))
								return npts;
							    
							hitp = 0;
						    } else {
							dot1 = dot(r_dir1, hits[i].hit_normal);
							if (dot1 >= 0.0)
							    continue;

							/* remember hits[i].hit_point as in point */
							hitp = &hits[i];
						    }
						}
					    }		

					}   

				else
				{
			    /* check for hits on the "foundation" pillar under the top.  The
			     * ray may have entered through the top of the pillar, possibly
			     * after having come down through the triangles above
			     */
			    bbmax.z = node->dsp_min[Z];
			    bbmin.z = 0.0;
			    if (dsp_in_rpp(r_pt1, invdir, bbmin, bbmax, &r_min, &r_max, &dmin, &dmax)) 
			    {
					/* hit rpp */
					hits[0].hit_vpriv = (double3) (0.0);
					hits[1].hit_vpriv = (double3) (0.0);

					minpt = r_pt1 + r_min * r_dir1;
					maxpt = r_pt1 + r_max * r_dir1;

					hits[0].hit_dist = r_min;
					hits[0].hit_surfno = dmin;
					hits[0].hit_point = vload3(0, minpt);
					hits[0].hit_normal = vload3(0, dsp_pl[dmin]);

					hits[1].hit_dist = r_max;
					hits[1].hit_surfno = dmax;
					hits[1].hit_point = vload3(0, maxpt);
					hits[1].hit_normal = vload3(0, dsp_pl[dmax]);

					/* add a segment to the list */
					
					npts+= dsp_segp(res, idx, &hits[0], &hits[1], dsp, r_dir);
					nhits++;
					if (( r_min > 0.0 || r_max> 0.0 )&&(nhits > a_onehit))
						return npts;
					    
			    }
				}
			}
	}

    return npts;
}






void compute_normal_at_gridpoint(double3 *N, const struct dsp_specific *dsp, int x, int y)
{
    double3 A, C, D, E, tmp;
    double3 Vac, Vde;

    tmp = double3(x, y, DSP(dsp, x, y));

    if (x == 0) {
	A = vload3(0, tmp);
    } else {
	A = double3(x-1, y, DSP(dsp, x-1, y));
    }

    if (x >= (dsp->dsp_xcnt-1)) {
	C = vload3(0, tmp);
    } else {
	C = double3(x+1, y,  DSP(dsp, x+1, y));
    }

    if (y == 0) {
	D = vload3(0, tmp);
    } else {
	D = double3(x, y-1, DSP(dsp, x, y-1));
    }

    if (y >= (dsp->dsp_ycnt-1)) {
	E = vload3(0, tmp);
    } else {
	E = double3(x, y+1, DSP(dsp, x, y+1));
    }


    /* Computing in world coordinates */
    tmp = vload3(0, A);
    A = MAT4X3PNT(vload16(0, dsp->dsp_stom), tmp);

    tmp = vload3(0, C);
    C = MAT4X3PNT(vload16(0, dsp->dsp_stom), tmp);

    tmp = vload3(0, D);
    D = MAT4X3PNT(vload16(0, dsp->dsp_stom), tmp);

    tmp = vload3(0, E);
    E = MAT4X3PNT(vload16(0, dsp->dsp_stom), tmp);

    Vac = C - vload3(0, A);
    Vde = E - vload3(0, D);

    Vac = normalize(Vac);
    Vde = normalize(Vde);
    *N = cross(Vac, Vde);

    *N = normalize(*N);
}




/**
 * Given ONE ray distance, return the normal and entry/exit point.
 */
void dsp_norm(struct hit *hitp, const double3 r_pt, const double3 r_dir, global const struct dsp_specific *dsp)
{
    double3 N, t, tmp, A;
    double3 Anorm, Bnorm, Dnorm, Cnorm, ABnorm, CDnorm;
    double Xfrac, Yfrac;
    int x, y;
    double3 pt;
    double dot1;
    double len;

    /* compute hit point */
    hitp->hit_point = r_pt + hitp->hit_dist * r_dir;


    if (hitp->hit_surfno != ZTOP || !dsp->dsp_smooth) 
    {
		/* We've hit one of the sides or bottom, or the user didn't
		 * ask for smoothing of the elevation data, so there's no
		 * interpolation to do.  Just transform the normal to model
		 * space, and compute the actual hit point
		 */

		/* transform normal into model space */
		tmp = MAT4X3VEC(dsp->dsp_mtos, hitp->hit_normal);
		tmp = normalize(tmp);
		hitp->hit_normal = vload3(tmp);
		
		return;
	}

	/* compute the distance between grid points in model space */
    tmp = double3(1.0, 0.0, 0.0);
    t = MAT4X3VEC(dsp->dsp_stom, tmp);
    len = length(t);

    /* get the cell we hit */
    x = hitp->hit_vpriv.x;
    y = hitp->hit_vpriv.y;

    compute_normal_at_gridpoint(&Anorm, dsp, x, y);
    compute_normal_at_gridpoint(&Bnorm, dsp, x+1, y);
    compute_normal_at_gridpoint(&Dnorm, dsp, x+1, y+1);
    compute_normal_at_gridpoint(&Cnorm, dsp, x, y+1);

    /* transform the hit point into DSP space for determining
     * interpolation
     */
    pt = MAT4X3PNT( vload16(0, dsp->dsp_mtos), hitp->hit_point);

    Xfrac = (pt.x - x);
    Yfrac = (pt.y - y);

    if (Xfrac < 0.0) Xfrac = 0.0;
    else if (Xfrac > 1.0) Xfrac = 1.0;

    if (Yfrac < 0.0) Yfrac = 0.0;
    else if (Yfrac > 1.0) Yfrac = 1.0;


    if (dsp->dsp_smooth == 2) {
	/* This is an experiment to "flatten" the curvature of the dsp
	 * near the grid points
	 */
#define SMOOTHSTEP(x)  ((x)*(x)*(3 - 2*(x)))
	Xfrac = SMOOTHSTEP(Xfrac);
	Yfrac = SMOOTHSTEP(Yfrac);
#undef SMOOTHSTEP
    }

    /* we compute the normal along the "X edges" of the cell */
    Anorm = Anorm * (1.0-Xfrac);
    Bnorm = Bnorm * Xfrac;
    ABnorm = Anorm + Bnorm;
    ABnorm = normalize(ABnorm);

    Cnorm = Cnorm * (1.0-Xfrac);
    Dnorm = Dnorm * Xfrac;
    CDnorm = Dnorm + Cnorm;
    CDnorm = normalize(CDnorm);

    /* now we interpolate the two X edge normals to get the final one
     */
    ABnorm = ABnorm * (1.0-Yfrac);
    CDnorm = CDnorm * Yfrac;
    N = ABnorm + CDnorm;

    N = normalize(N);

    dot1 = dot(N, r_dir);

    if ((ZERO(hitp->hit_vpriv.z) && dot1 > 0.0)/* in-hit needs fix */ ||
	(ZERO(hitp->hit_vpriv.z - 1.0) && dot1 < 0.0)/* out-hit needs fix */) {
	/* bring the normal back to being perpendicular to the ray to
	 * avoid "flipped normal" warnings
	 */
	A = cross(r_dir, N);
	N = cross(A, r_dir);
	N = normalize(N);
	dot1 = dot(N, r_dir);
    }

    hitp->hit_normal = vload3(N);
}

#undef XMIN 
#undef XMAX 
#undef YMIN 
#undef YMAX 
#undef ZMIN 
#undef ZMAX 
#undef ZMID 
#undef ZTOP 
#undef DSP_CUT_DIR_ADAPT
#undef DSP_CUT_DIR_llUR 
#undef DSP_CUT_DIR_ULlr 



/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */

