#include "common.cl"

#define SQ(a) ((a)*(a))
#define METABALL_METABALL 0
#define METABALL_ISOPOTENTIAL 1
#define METABALL_BLOB 2

struct wdb_metaballpt {
    double fldstr; /**< @brief  field strength */
    double sweat;  /**< @brief  beta value used for metaball and blob evaluation */
    double3 coord;
};

struct metaball_specific {
    /* these three defines are used with the method field */
    int method;
    double threshold;
    double initstep;
    double finalstep; /* for raytrace stepping. */
    struct wdb_metaballpt *metaball_ctrl_head;
};

double rt_metaball_point_value_metaball(const double3 *p, const struct wdb_metaballpt *points)
{
    /* Makes the compiler happy */
    return 0.0;
}


double
rt_metaball_point_value_iso(const double3 *p, const struct wdb_metaballpt *points)
{
    struct wdb_metaballpt *mbpt;
    double ret = 0.0;
    double3 v;

    for(*mbpt = *points; *mbpt!=NULL; mbpt++) {
	v = mbpt->coord - vload3(0, *p);
	ret += fabs(mbpt->fldstr) * mbpt->fldstr / dot(v, v);	/* f/r^2 */
	}
    return ret;
}


double
rt_metaball_point_value_blob(const double3 *p, const struct wdb_metaballpt *points)
{
    struct wdb_metaballpt *mbpt;
    double ret = 0.0;
    double3 v;

    for(*mbpt = *points; *mbpt!=NULL; mbpt++)  {
		/* TODO: test if sweat is sufficient enough that r=0 returns a positive value? */
		/* TODO: test to see if negative contribution needs to be wiped out? */
		v = mbpt->coord - vload3(0, *p);
		ret += 1.0 / exp((mbpt->sweat/(mbpt->fldstr*mbpt->fldstr)) * dot(v, v) - mbpt->sweat);
    }
    return ret;
}


/* main point evaluation function, to be exposed to the ugly outside
 * world.
 */
double
rt_metaball_point_value(const double3 *p, const struct metaball_specific *mb)
{
    switch (mb->method) {
	case METABALL_METABALL:
	    return rt_metaball_point_value_metaball(p, &mb->metaball_ctrl_head);
	case METABALL_ISOPOTENTIAL:
	    return rt_metaball_point_value_iso(p, &mb->metaball_ctrl_head);
	case METABALL_BLOB:
	    return rt_metaball_point_value_blob(p, &mb->metaball_ctrl_head);
	default:
	    break;
    }
    return 0;
}

/*
 * Solve the surface intersection of mb with an accuracy of finalstep given that
 * one of the two points (a and b) are inside and the other is outside.
 */
int
rt_metaball_find_intersection(double3 *intersect, const struct metaball_specific *mb, double3 *a, double3 *b, double step, const double finalstep)
{
    double3 mid;
    const double3 *midp = (const double3 *)&mid;

    do
    {	
    	mid = *a + *b;
    	mid = mid * 0.5;
    	if ((rt_metaball_point_value(a, mb) >= mb->threshold) != (rt_metaball_point_value(midp, mb) >= mb->threshold))
    		{ b = a; }
    	*a = mid;
    	step = step/2.0;
    }while(finalstep<=step);

	*intersect = vload3(0, mid);	/* should this be the midpoint between a and b? */
	return 0;
}

int 
metaball_shot(RESULT_TYPE *res, const double3 r_pt, const double3 r_dir, const double r_min, const double r_max, const uint idx, global const struct metaball_specific *mb)
{
    struct hit hits[2];
    int retval = 0;
    int a_onehit = 1;
    double step, distleft;
    double3 p, inc;
    const double3 *cp = (const double3 *)&p;

    step = mb->initstep;
    distleft = (r_max-r_min) + step * 3.0;

    p = vload3(0, r_pt);
    inc = r_dir * step; /* assume it's normalized and we want to creep at step */

    /* walk back out of the solid */
    while (rt_metaball_point_value(cp, mb) >= mb->threshold) {
	distleft += step;
	p = p - vload3(0, inc);
    }

    
	int mb_stat = 0, segsleft = fabs(a_onehit);
	double3 lastpoint;

	while (distleft >= 0.0 || mb_stat == 1) {
	    /* advance to the next point */
	    distleft -= step;
	    lastpoint = vload3(0, p);
	    p = p + vload3(0, inc);
	    if (mb_stat == 1) {
		if (rt_metaball_point_value(cp, mb) < mb->threshold) {
		    double3 intersect, delta;
		    const double3 *pA = (const double3 *)&lastpoint;
		    const double3 *pB = (const double3 *)&p;
		    rt_metaball_find_intersection(&intersect, mb, pA, pB, step, mb->finalstep);
		    hits[1].hit_point = vload3(0, intersect);
		    --segsleft;
		    ++retval;
		    delta = intersect -vload3(0, r_pt);
		    hits[1].hit_dist = length(delta);
		    hits[1].hit_surfno = 0;
		    mb_stat = 0;
		    if (a_onehit != 0 && segsleft <= 0)
				return retval;
		}
	    } else {
		if (rt_metaball_point_value(cp, mb) > mb->threshold) {
		    double3 intersect, delta;
		    const double3 *pA = (const double3 *)&lastpoint;
		    const double3 *pB = (const double3 *)&p;
		    rt_metaball_find_intersection(&intersect, mb, pA, pB, step, mb->finalstep);
		    --segsleft;
		    ++retval;
		    hits[0].hit_point = vload3(0, intersect);
		    delta = intersect - vload3(0, r_pt);
		    hits[0].hit_dist = length(delta);
		    hits[0].hit_surfno = 0;
		    do_segp(res, idx, &hits[0], &hits[1]);

		    mb_stat = 1;
		    step = mb->initstep;
		}
	    }
	}
    return retval;
}


/**
 * Given ONE ray distance, return the normal and entry/exit point.
 */
void metaball_norm(struct hit *hitp, const double3 r_pt, const double3 r_dir, global const struct metaball_specific *mb)
{
    struct wdb_metaballpt *mbpt;
    double3 v;
    double a;

    hitp->hit_normal = (double3)(0.0);

    switch (mb->method) {
	case METABALL_METABALL:
	    break;
	case METABALL_ISOPOTENTIAL:
	    for(mbpt = mb->metaball_ctrl_head; *mbpt!=NULL; mbpt++)  {
		v = hitp->hit_point - vload3(0, mbpt->coord);
		a = dot(v, v);
		hitp->hit_normal = hitp->hit_normal + (fabs(mbpt->fldstr)*mbpt->fldstr / (a*a)) * v;	/* f/r^4 */
	    }
	    break;
	case METABALL_BLOB:
	    for(mbpt = mb->metaball_ctrl_head; *mbpt!=NULL; mbpt++)  {
		v = hitp->hit_point - vload3(0, mbpt->coord);
		a = dot(v, v);
		hitp->hit_normal = hitp->hit_normal + 2.0*mbpt->sweat/SQ(mbpt->fldstr)*exp(mbpt->sweat*(1-(a/SQ(mbpt->fldstr)))) * v;
	    }
	    break;
	default: break;
    }
    hitp->hit_normal = normalize(hitp->hit_normal);

    return;
}
