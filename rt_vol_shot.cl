#include "common.cl"

struct rt_vol_specific {
    double vol_xnorm[3];	/* local +X norm in model coords */
    double vol_ynorm[3];
    double vol_znorm[3];
    double vol_mat[16];	/* model to ideal space */
    double vol_origin[3];	/* local coords of grid origin (0, 0, 0) for now */
    double vol_large[3];	/* local coords of XYZ max */
    double vol_cellsize[3];
    uint xdim;
    uint ydim;
    uint zdim;
    uint lo;        /**< @brief Low threshold */
    uint hi;        /**< @brief High threshold */
    uchar map[1];
};

/*
 * Codes to represent surface normals.
 * In a bitmap, there are only 4 possible normals.
 * With this code, reverse the sign to reverse the direction.
 * As always, the normal is expected to point outwards.
 */
#define NORM_ZPOS 3
#define NORM_YPOS 2
#define NORM_XPOS 1
#define NORM_XNEG (-1)
#define NORM_YNEG (-2)
#define NORM_ZNEG (-3)

#define X 0
#define Y 1
#define Z 2

/*
 * Regular bit addressing is used:  (0..W-1, 0..N-1),
 * but the bitmap is stored with two cells of zeros all around,
 * so permissible subscripts run (-2..W+1, -2..N+1).
 * This eliminates special-case code for the boundary conditions.
 */
#define VOL_XWIDEN 2
#define VOL_YWIDEN 2
#define VOL_ZWIDEN 2
#define VOL(_vip, _xx, _yy, _zz)	(vip)->map[]

static inline uint VOL(global const struct vol_specific *vip, uint x, uint y, uint z) {
    uint id = ((z + VOL_ZWIDEN) * (vip->ydim + VOL_YWIDEN*2)+  (y + VOL_YWIDEN)) * (vip->xdim + VOL_XWIDEN*2) + x + VOL_XWIDEN ;
    return ((vip->map[(id >> 3)] & (1 << (id & 7))));
}

/**
 * Transform the ray into local coordinates of the volume ("ideal space").
 * Step through the 3-D array, in local coordinates.
 * Return intersection segments.
 *
 */
int rt_vol_shot(RESULT_TYPE *res, const double3 r_pt_, const double3 r_dir_, const uint idx, global const struct rt_vol_specific *volp)
{
    struct hit hits[2];
    double3 invdir;
    double t0;	/* in point of cell */
    double t1;	/* out point of cell */
    double tmax;	/* out point of entire grid */
    double t[3];	/* next t value for XYZ cell plane intersect */
    double delta[3];	/* spacing of XYZ cell planes along ray */
    int igrid[3];/* Grid cell coordinates of cell (integerized) */
    double3 P;	/* hit point */
    int inside;	/* inside/outside a solid flag */
    int in_axis;
    int out_axis;
    int j;
    int npts, nhits;
    double3 r_pt, dir;
    double r_dir[3];
    int rt_vol_normtab[3] = { NORM_XPOS, NORM_YPOS, NORM_ZPOS };

    /* Transform actual ray into ideal space at origin in X-Y plane */
    r_pt = MAT4X3PNT(vload16(0, volp->vol_mat), r_pt_);
    dir = MAT4X3VEC(vload16(0, volp->vol_mat), r_dir_);
    
    /* compute the inverse of the direction cosines */
    if (!ZERO(dir.x)) {
	invdir.x = 1.0/dir.x;
	r_dir[X] = dir.x;
    } else {
	invdir.x = INFINITY;
	r_dir[X] = 0.0;
    }
    if (!ZERO(dir.y)) {
	invdir.y = 1.0/dir.y;
	r_dir[Y] = dir.y;
    } else {
	invdir.y = INFINITY;
	r_dir[Y] = 0.0;
    }
    if (!ZERO(dir.z)) {
	invdir.z = 1.0/dir.z;
	r_dir[Z] = dir.z;
    } else {
	invdir.z = INFINITY;
	r_dir[Z] = 0.0;
    }

    /* intersect ray with ideal grid rpp */
    if (! rt_in_rpp2(r_pt, invdir, volp->vol_origin, volp->vol_large, &t0, &tmax))
	    return 0;	/* MISS */
    
    P = r_pt + t0 * vload3(0, r_dir);	/* P is hit point */

    /* find grid cell where ray first hits ideal space bounding RPP */
    igrid[X] = (P.x - volp->vol_origin[X]) / volp->vol_cellsize[X];
    igrid[Y] = (P.y - volp->vol_origin[Y]) / volp->vol_cellsize[Y];
    igrid[Z] = (P.z - volp->vol_origin[Z]) / volp->vol_cellsize[Z];
    if (igrid[X] < 0)
	    igrid[X] = 0;
    else if (igrid[X] >= volp->xdim) 
	    igrid[X] = volp->xdim-1;
    
    if (igrid[Y] < 0) 
	    igrid[Y] = 0;
    else if (igrid[Y] >= volp->ydim) 
	    igrid[Y] = volp->ydim-1;
    
    if (igrid[Z] < 0) 
	    igrid[Z] = 0;
    else if (igrid[Z] >= volp->zdim) 
	    igrid[Z] = volp->zdim-1;
    

    /* X setup */
    if (ZERO(r_dir[X])) {
	t[X] = INFINITY;
	delta[X] = 0;
    } else {
	j = igrid[X];
	if (r_dir[X] < 0) j++;
	t[X] = (volp->vol_origin[X] + j*volp->vol_cellsize[X] - r_pt.x) * invdir.x;
	delta[X] = volp->vol_cellsize[X] * fabs(invdir.x);
    }
    /* Y setup */
    if (ZERO(r_dir[Y])) {
	t[Y] = INFINITY;
	delta[Y] = 0;
    } else {
	j = igrid[Y];
	if (r_dir[Y] < 0) j++;
	t[Y] = (volp->vol_origin[Y] + j*volp->vol_cellsize[Y] - r_pt.y) * invdir.y;
	delta[Y] = volp->vol_cellsize[Y] * fabs(invdir.y);
    }
    /* Z setup */
    if (ZERO(r_dir[Z])) {
	t[Z] = INFINITY;
	delta[Z] = 0;
    } else {
	j = igrid[Z];
	if (r_dir[Z] < 0) j++;
	t[Z] = (volp->vol_origin[Z] + j*volp->vol_cellsize[Z] - r_pt.z) * invdir.z;
	delta[Z] = volp->vol_cellsize[Z] * fabs(invdir.z);
    }

    /* Find face of entry into first cell -- max initial t value */
    if (t[X] >= t[Y]) {
	in_axis = X;
	t0 = t[X];
    } else {
	in_axis = Y;
	t0 = t[Y];
    }
    if (t[Z] > t0) {
	in_axis = Z;
	t0 = t[Z];
    }
    
    /* Advance to next exits */
    t[X] += delta[X];
    t[Y] += delta[Y];
    t[Z] += delta[Z];

    /* Ensure that next exit is after first entrance */
    if (t[X] < t0) 
	    t[X] += delta[X];

    if (t[Y] < t0) 
	    t[Y] += delta[Y];
    
    if (t[Z] < t0) 
	    t[Z] += delta[Z];
    
    npts = 0;
    nhits = 0;
    while (t0 < tmax) {
	int val;

	/* find minimum exit t value */
	if (t[X] < t[Y]) {
	    if (t[Z] < t[X]) {
		out_axis = Z;
		t1 = t[Z];
	    } else {
		out_axis = X;
		t1 = t[X];
	    }
	} else {
	    if (t[Z] < t[Y]) {
		out_axis = Z;
		t1 = t[Z];
	    } else {
		out_axis = Y;
		t1 = t[Y];
	    }
	}

	/* Ray passes through cell igrid[XY] from t0 to t1 */
	val = VOL(volp, igrid[X], igrid[Y], igrid[Z]);

	if (!(nhits&1)) {
	    if ((val >= volp->lo) && (val <= volp->hi)) {
		/* Handle the transition from vacuum to solid */
		/* Start of segment (entering a full voxel) */
		hits[0].hit_dist = t0;
		/* Compute entry normal */
		hits[0].hit_surfno = (r_dir[in_axis] < 0) ? rt_vol_normtab[in_axis] : (-rt_vol_normtab[in_axis]);
		nhits++;
	    } else {
            /* Do nothing, marching through solid */
        }
	} else {
	    if ((val >= volp->lo) && (val <= volp->hi)) {
		/* Do nothing, marching through solid */
	    } else {
		/* End of segment (now in an empty voxel) */
		/* Handle transition from solid to vacuum */
		hits[1].hit_dist = t0;
		/* Compute exit normal */
		hits[1].hit_surfno = (r_dir[in_axis] < 0) ?	(-rt_vol_normtab[in_axis]) : rt_vol_normtab[in_axis];
		nhits++;
        npts += 2;
        do_segp(res, idx, &hits[0], &hits[1]);
	    }
	}

	/* Take next step */
	t0 = t1;
	in_axis = out_axis;
	t[out_axis] += delta[out_axis];
	if (r_dir[out_axis] > 0)
	    igrid[out_axis]++;
	else 
	    igrid[out_axis]--;
    }

    if ((nhits&1)) {
		/* close off the final segment */
		hits[1].hit_dist = tmax;
        /* Compute exit normal.  Previous out_axis is now in_axis */
        hits[1].hit_surfno = (r_dir[in_axis] < 0) ? (-rt_vol_normtab[in_axis]) : rt_vol_normtab[in_axis];
		nhits++;
        npts += 2;
        do_segp(res, idx, &hits[0], &hits[1]);
	}
    return npts;
}


/**
 * Given one ray distance, return the normal and
 * entry/exit point.
 * This is mostly a matter of translating the stored
 * code into the proper normal.
 */
void rt_vol_norm(struct hit *hitp, const double3 r_pt, const double3 r_dir, global const struct rt_vol_specific *volp)
{
    hitp->hit_point = r_pt + r_dir * hitp->hit_dist;

    switch (hitp->hit_surfno) {
	case NORM_XPOS: hitp->hit_normal =  vload3(0, volp->vol_xnorm); break;
	case NORM_XNEG: hitp->hit_normal = -vload3(0, volp->vol_xnorm); break;
	case NORM_YPOS: hitp->hit_normal =  vload3(0, volp->vol_ynorm); break;
	case NORM_YNEG: hitp->hit_normal = -vload3(0, volp->vol_ynorm); break;
	case NORM_ZPOS: hitp->hit_normal =  vload3(0, volp->vol_znorm); break;
	case NORM_ZNEG: hitp->hit_normal = -vload3(0, volp->vol_znorm); break;
	default: hitp->hit_normal = (double3)(0.0); break;
    }
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