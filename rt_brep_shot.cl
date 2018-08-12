#include "common.cl"

#define BREP_GRAZING_DOT_TOL 0.000017453

#define	CLEAN_HIT 1
#define	NEAR_HIT 2
#define	NEAR_MISS 3
#define	CRACK_HIT 4 //applied to first point of two near_miss points with same normal direction, second point removed
    
#define	ENTERING 1
#define	LEAVING -1

#define MAX_BREP_SUBDIVISION_INTERSECTS 5
#define BREP_INTERSECT_ROOT_DIVERGED -2
#define BREP_INTERSECT_FOUND 1

#define BREP_INTERSECTION_ROOT_EPSILON 1e-6
#define BREP_MAX_ITERATIONS 100

#define ON_UNSET_VALUE -1.23432101234321e+308
#define ON_DBL_MIN 2.22507385850720200e-308
#define DBL_MAX 1.7976931348623157E+308
#define ON_EPSILON 2.2204460492503131e-16
#define ON_SQRT_EPSILON 1.490116119385000000e-8

#define SMALL_FASTF 1.0e-77
#define ROOT_TOL 1.0e-7
#define BREP_EDGE_MISS_TOLERANCE 5e-3
#define VUNITIZE_TOL 1.0e-15


struct Stl{
	struct BBNode* m_children;
	struct BRNode* m_trims_above;
};


struct BBNode{	
	double m_min[3];
	double m_max[3];
    double m_u[2];
    double m_v[2];
    bool m_checkTrim;
    bool m_trimmed;
    double3 m_normal;
    struct Stl *m_stl;
    struct ON_BrepFace m_face;
};


struct BRNode{
	double m_min[3];
	double m_max[3];
    int m_adj_face_index;
    bool m_XIncreasing;
    bool m_Vertical;
    double m_u[2];
    double m_v[2];
    double m_t[2];
    bool m_trimmed;
    double3 m_start;
    double3 m_end;
    struct ON_Curve m_trim;
};



struct brep_hit{
    struct ON_BrepFace face;
    double dist;
    double3 origin;
    double3 point;
    double3 normal;
    double3 uv;
    int trimmed;
    bool closeToEdge;
    int oob;
    int hit;
    int direction;
    int m_adj_face_index;
    // XXX - calculate the dot of the dir with the normal here!
    struct BBNode *sbv;
    int useful;
};


struct brep_specific
{
	struct BBNode bvh;
};


void Surface_Ev1Der( // returns false if unable to evaluate
       double s, double t, // evaluation parameters
       double3* point, double3* ds, double3* dt,
       int side,        // optional - determines which side to evaluate from
                       //         0 = default
                       //         1 from NE quadrant
                       //         2 from NW quadrant
                       //         3 from SW quadrant
                       //         4 from SE quadrant
       int* hint       // optional - evaluation hint (int[2]) used to speed
                       //            repeated evaluations
       )
{
  bool rc = false;
  int dim = Dimension();
  double ws[3*32];
  double* v;

  point = (double3)(0.0);
  *ds = (double3) (0.0);
  *dt = (double3)(0.0);

  if ( dim <= 32 )
    v = ws;
  
  else 
    v = (double*)onmalloc(3*dim*sizeof(*v));
  
  rc = Evaluate( s, t, 1, dim, v, side, hint );
  *point.x = v[0];
  *ds.x = v[dim];
  *dt.x = v[2*dim];
  if ( dim > 1 ) {
    *point.y = v[1];
    *ds.y = v[dim+1];
    *dt.y = v[2*dim+1];
    if ( dim > 2 ) {
      *point.z = v[2];
      *ds.z = v[dim+2];
      *dt.z = v[2*dim+2];
      if ( dim > 32 )
        onfree(v);
    }
  }
  return rc;
}



bool Curve_Ev1Der( // returns false if unable to evaluate
       double t,         // evaluation parameter
       double3* point, double3* derivative,
       int side,        // optional - determines which side to evaluate from
                       //      <= 0 to evaluate from below, 
                       //      >  0 to evaluate from above
       int* hint       // optional - evaluation hint used to speed repeated
                       //            evaluations
       )
{
  bool rc = false;
  int dim = Dimension();
  double ws[2*64];
  double* v;
  *point = (double3) (0.0);
  *derivative = (double3) (0.0);

  if ( dim <= 64 ) 
    v = ws;
  
  else
    v = (double*)onmalloc(2*dim*sizeof(*v));
  
  rc = Evaluate( t, 1, dim, v, side, hint );
  *point.x = v[0];
  *derivative.x = v[dim];
  if ( dim > 1 ) {
    *point.y = v[1];
    *derivative.y = v[dim+1];
    if ( dim > 2 ) {
      *point.z = v[2];
      *derivative.z = v[dim+2];
      if ( dim > 64 )
        onfree(v);
    }
  }
  return rc;
}


bool Ev2Der( // returns false if unable to evaluate
       double t,         // evaluation parameter
       double3* point, double3* firstDervative, double3* secondDervative,
       int side,        // optional - determines which side to evaluate from
                       //      <= 0 to evaluate from below, 
                       //      >  0 to evaluate from above
       int* hint       // optional - evaluation hint used to speed repeated
                       //            evaluations
       )
{
  bool rc = false;
  int dim = Dimension();
  double ws[3*64];
  double* v;
  *point = (double3) (0.0);
  *firstDervative = (double3) (0.0);
  *secondDervative = (double3)(0.0);
  
  if ( dim <= 64 )
    v = ws;
  
  else 
    v = (double*)onmalloc(3*dim*sizeof(*v));
  
  rc = Evaluate( t, 2, dim, v, side, hint );
  *point.x = v[0];
  *firstDervative.x = v[dim];
  *secondDervative.x = v[2*dim];
  if ( dim > 1 ) {
    *point.y = v[1];
    *firstDervative.y = v[dim+1];
    *secondDervative.y = v[2*dim+1];
    if ( dim > 2 ) {
      *point.z = v[2];
      *firstDervative.z = v[dim+2];
      *secondDervative.z = v[2*dim+2];
      if ( dim > 64 )
        onfree(v);
    }
  }

  return rc;
}


double ParameterAt(double* m_t, double x)  
{
  return ((x != ON_UNSET_VALUE && isfinite(x)) ? ((1.0-x)*m_t[0] + x*m_t[1]) : ON_UNSET_VALUE);
}


bool IsIncreasing(double* m_t)
{
  return ( (m_t[0] < m_t[1]) && (m_t[0] != ON_UNSET_VALUE) && (isfinite(m_t[0])) && (m_t[1] != ON_UNSET_VALUE) && (isfinite(m_t[1])) )  ? true : false;
}


// returns tminus < tplus: parameters tminus <= s <= tplus
bool GetParameterTolerance(double t, double* tminus, double* tplus)
{
  bool rc = false;
  ON_Interval d = Domain();
  double dt;
  if (IsIncreasing(d))
  {
  rc = (d[0] < d[1]) ? true : false;
	  if ( rc ) {
	    if ( t < d[0] )
	      t = d[0];
	    else if (t > d[1] )
	      t = d[1];
	    dt = (d[1]-d[0])*8.0* ON_SQRT_EPSILON + (fabs(d[0]) + fabs(d[1]))* ON_EPSILON;
	    if (dt >= (d[1]-d[0]))
	      dt = 0.5*(d[1]-d[0]);
	    if ( tminus )
	      *tminus = t-dt;
	    if ( tplus )
	      *tplus = t+dt;
	 	 }
	return rc;
	}   
}



bool EvTangent(double t, double3* point, double3* tangent, int side, int* hint)
{
  double3 D1, D2;
  *tangent = (double3)(0.0);
  int rc = Curve_Ev1Der(t, point, tangent, side, hint);
  *tangent = normalize(*tangent);
  if (rc && !(*tangent.x == 0 && *tangent.y==0 && *tangent.z==0)) 
  {
    if (Ev2Der(t, point, &D1, &D2, side, hint))
    {
      // Use l'Hopital's rule to show that if the unit tanget
      // exists, the 1rst derivative is zero, and the 2nd
      // derivative is nonzero, then the unit tangent is equal
      // to +/-the unitized 2nd derivative.  The sign is equal
      // to the sign of D1(s) o D2(s) as s approaches the 
      // evaluation parameter.
      *tangent = normalize(D2);
      rc = (*tangent.x == 0 && *tangent.y==0 && *tangent.z==0) ?0 :1;
      if (rc)
      {
        ON_Interval domain = Domain();
        double tminus = 0.0;
        double tplus = 0.0;
        if (IsIncreasing(domain) && GetParameterTolerance(t, &tminus, &tplus))
        {
          double3 p;
          double3 d1, d2;
          double eps = 0.0;
          double d1od2tol = 0.0; //1.0e-10; // 1e-5 is too big
          double d1od2;
          double tt = t;

          if ( (t < domain[1] && side >= 0) || (t == domain[0]) )
          {
            eps = tplus-t;
            if ( eps <= 0.0 || t+eps > ParameterAt(domain, 0.1) )
              return rc;
          }
          else if ( (t > domain[0] && side < 0) || (t == domain[1]) )
          {
            eps = tminus - t;
            if ( eps >= 0.0 || t+eps < ParameterAt(domain, 0.9) )
              return rc;
          }

          int i, negative_count=0, zero_count=0;
          int test_count = 3;
          for ( i = 0; i < test_count; i++, eps *= 0.5 )
          {
            tt = t + eps;
            if (tt==t)
              break;
            if (!Ev2Der(tt, &p, &d1, &d2, side, 0))
              break;
            d1od2 = d1*d2;
            if (d1od2 > d1od2tol)
              break;
            if (d1od2 < d1od2tol)
              negative_count++;
            else
              zero_count++;
          }
          if ( negative_count > 0 && test_count == negative_count+zero_count )
          		*tangent = -(*tangent);
          
        }
      }
    }
  }
  return rc;
}


bool EvPoint( // returns false if unable to evaluate
       double t,         // evaluation parameter
       double3* point,   // returns value of curve
       ON_Curve m_trim,
       int side,        // optional - determines which side to evaluate from
                       //         0 = default
                       //      <  0 to evaluate from below, 
                       //      >  0 to evaluate from above
       int* hint       // optional - evaluation hint used to speed repeated
                       //            evaluations
       ) 
{
  bool rc = false;
  double ws[128];
  double* v;
  if (Dimension() <= 3) {
    v = *point.x;
    *point = (double3)(0.0);
  }
  else if ( Dimension() <= 128 )
    v = ws;
 
  else 
    v = (double*)onmalloc(Dimension()*sizeof(*v));
 
  rc = Evaluate(t, 0, Dimension(), v, side, hint);
  if (Dimension() > 3) {
    *point.x = v[0];
    *point.y = v[1];
    *point.z = v[2];
    if(Dimension() > 128)
      onfree(v);
  }
  return rc;
}


double BRNode_getCurveEstimateOfV(struct BRNode br, double u, double tol)
{
    double3 tangent;
    double3 A, B;
    double Ta, Tb;
    double3 Tan_start, Tan_end;
    double3 point;

    br.m_start = (double3)(0.0);
    if (!EvPoint(br.m_t[0], br.m_start, br.m_trim, 0, 0))
    	br.m_start = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

    br.m_end = (double3)(0.0);
    if (!EvPoint(br.m_t[1], br.m_end, br.m_trim, 0, 0))
    	br.m_end = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

    if (br.m_start.x < br.m_end.x) {
	A = vload3(0, br.m_start);
	B = vload3(0, br.m_end);
	Ta = br.m_t[0];
	Tb = br.m_t[1];
    } else {
	A = vload3(0, br.m_end);
	B = vload3(0, br.m_start);
	Ta = br.m_t[1];
	Tb = br.m_t[0];
    }

    double dU = B.x - A.x;
    if (NEAR_ZERO(dU, tol))    /* vertical */
		return A.y;
    
    EvTangent(Ta, &point, &Tan_start, 0, 0);
    EvTangent(Tb, &point, &Tan_end, 0, 0);

    double dT = Tb - Ta;
    double guess;
    double3 p = (double3)(0.0);

    /* Use quick binary subdivision until derivatives at end points in 'u' are within 5 percent */
    while (!NEAR_ZERO(dU, tol) && !NEAR_ZERO(dT, tol)) {
	guess = Ta + dT / 2;
	if ( !EvPoint(guess, p, br.m_trim, 0, 0) )
    	p = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

	if (UNLIKELY(NEAR_EQUAL(p.x, u, SMALL_FASTF))) /* hit 'u' exactly, done deal */
	    return p.y;

	if (p.x > u) {
	    /* v is behind us, back up the end */
	    Tb = guess;
	    B = vload3(0, p);
	    EvTangent(Tb, &point, &Tan_end, 0, 0);
	} else {
	    /* v is in front, move start forward */
	    Ta = guess;
	    A = vload3(0, p);
	    EvTangent(Ta, &point, &Tan_start, 0, 0);
	}
	dT = Tb - Ta;
	dU = B.x - A.x;
    }

    dU = B.x - A.x;
    if (NEAR_ZERO(dU, tol))    /* vertical */
		return A.y;

    guess = Ta + (u - A.x) * dT / dU;
	if ( !EvPoint(guess, p, br.m_trim, 0, 0) )
    	p = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

    int cnt = 0;
    while ((cnt < 1000) && (!NEAR_EQUAL(p.x, u, tol))) {
	if (p.x < u) {
	    Ta = guess;
	    A = vload3(0, p);
	} else {
	    Tb = guess;
	    B = vload3(0, p);
	}
	dU = B.x - A.x;
	if (NEAR_ZERO(dU, tol))  /* vertical */
	    return A.y;

	dT = Tb - Ta;
	guess = Ta + (u - A.x) * dT / dU;
	if ( !EvPoint(guess, p, br.m_trim, 0, 0) )
    	p = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);
	cnt++;
    }
    return p.y;
}





bool BRNode_isTrimmed(struct BRNode br, double3 uv, double *trimdist)
{
    double3 bmin, bmax;
	double3 UnsetPoint = (double3) (ON_UNSET_VALUE,ON_UNSET_VALUE,ON_UNSET_VALUE);
 
    bmin.x = INFINITY;
    bmin.y = INFINITY;
    bmin.z = INFINITY;
	
	bmax.x = -INFINITY;
	bmax.y = -INFINITY;
	bmax.z = -INFINITY;

	br.m_start = (double3)(0.0);
    if ( !EvPoint(br.m_t[0], br.m_start, br.m_trim, 0, 0) )
    	br.m_start = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

    br.m_end = (double3)(0.0);
    if ( !EvPoint(br.m_t[1], br.m_end, br.m_trim, 0, 0) )
    	br.m_end = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

	if (br.m_start != UnsetPoint) 
	{
		bmin.x = fmin(bmin.x, br.m_start.x);
		bmin.y = fmin(bmin.y, br.m_start.y);
		bmin.z = fmin(bmin.z, br.m_start.z);

		bmax.x = fmax(bmax.x, br.m_start.x);
		bmax.y = fmax(bmax.y, br.m_start.y);
		bmax.z = fmax(bmax.z, br.m_start.z);
	}
	if (br.m_end != UnsetPoint) 
	{
		bmin.x = fmin(bmin.x, br.m_end.x);
		bmin.y = fmin(bmin.y, br.m_end.y);
		bmin.z = fmin(bmin.z, br.m_end.z);

		bmax.x = fmax(bmax.x, br.m_end.x);
		bmax.y = fmax(bmax.y, br.m_end.y);
		bmax.z = fmax(bmax.z, br.m_end.z);
	}

    if ((bmin.x <= uv.x) && (uv.x <= bmax.x))   /* if check trim and in BBox */
    {
		double v = BRNode_getCurveEstimateOfV(br, uv.x, 0.0000001);
		*trimdist = v - uv.y;
		if (uv.y <= v) {
		    if (br.m_XIncreasing)
				return true;
		    else
				return false;
		} else if (uv.y > v) {
		    if (!br.m_XIncreasing)
				return true;
		    else 
				return false;
		} else 
		    return true;
    } else {
	*trimdist = -1.0;
	if (br.m_trimmed) 
	    return true;
	else 
	    return false;
    }
}


int BBNode_getTrimsAbove(struct BBNode sbv, double3 uv, struct BRNode *out_leaves)
{
    double3 bmin, bmax;
    double dist;
    int len = 0;
    double3 UnsetPoint = (double3) (ON_UNSET_VALUE,ON_UNSET_VALUE,ON_UNSET_VALUE);

    for (int i = 0; i != m_stl->m_trims_above.end(); i++) {
	struct BRNode br = sbv.m_stl->m_trims_above[i];

    bmin.x = INFINITY;
    bmin.y = INFINITY;
    bmin.z = INFINITY;
	
	bmax.x = -INFINITY;
	bmax.y = -INFINITY;
	bmax.z = -INFINITY;

	br.m_start = (double3)(0.0);
    if ( !EvPoint(br.m_t[0], br.m_start, br.m_trim, 0, 0) )
    	br.m_start = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

    br.m_end = (double3)(0.0);
    if ( !EvPoint(br.m_t[1], br.m_end, br.m_trim, 0, 0) )
    	br.m_end = (double3)(ON_UNSET_VALUE, ON_UNSET_VALUE, ON_UNSET_VALUE);

	if (br.m_start != UnsetPoint) 
	{
		bmin.x = fmin(bmin.x, br.m_start.x);
		bmin.y = fmin(bmin.y, br.m_start.y);
		bmin.z = fmin(bmin.z, br.m_start.z);

		bmax.x = fmax(bmax.x, br.m_start.x);
		bmax.y = fmax(bmax.y, br.m_start.y);
		bmax.z = fmax(bmax.z, br.m_start.z);
	}
	if (br.m_end != UnsetPoint) 
	{
		bmin.x = fmin(bmin.x, br.m_end.x);
		bmin.y = fmin(bmin.y, br.m_end.y);
		bmin.z = fmin(bmin.z, br.m_end.z);

		bmax.x = fmax(bmax.x, br.m_end.x);
		bmax.y = fmax(bmax.y, br.m_end.y);
		bmax.z = fmax(bmax.z, br.m_end.z);
	}

	dist = 0.000001; /* 0.03*DIST_PT_PT(bmin, bmax); */
	if ((uv.x > bmin.x - dist) && (uv.x < bmax.x + dist))
	    	out_leaves[len++] = br;
    }
    return len;
}



bool BBNode_isTrimmed(struct BBNode sbv, double3 uv, struct BRNode **closest, double* closesttrim, double within_distance_tol)
{
	int len;
	int i;
	double dist;
	double v;
    struct BRNode *list;
    struct BRNode *trims;

    *closesttrim = -1.0;
    if (sbv.m_checkTrim) {
	len = BBNode_getTrimsAbove(sbv, uv, trims);

	if (len==0)
	    return true;
	else { /* find closest BB */
	    struct BRNode *vclosest = NULL;
	    struct BRNode *uclosest = NULL;
	    double currHeight = 0.0;
	    bool currTrimStatus = false;
	    bool verticalTrim = false;
	    bool underTrim = false;
	    double vdist = 0.0;
	    double udist = 0.0;
	    list = &trims[0];

	    for (i = 0; i < len; i++) {
		/* skip if trim below */
		if (list[i].m_max[1] + within_distance_tol < uv.y) 
		    continue;
		
		if (list[i].m_Vertical) {
		    if ((list[i].m_v[0] <= uv.y) && (list[i].m_v[1] >= uv.y)) {
			dist = fabs(uv.x - list[i].m_v[0]);
			if (!verticalTrim) { /* haven't seen vertical trim yet */
			    verticalTrim = true;
			    vdist = dist;
			    *vclosest = list[i];
			} else {
			    if (dist < vdist) {
				vdist = dist;
				*vclosest = list[i];
			    }
			}

		    }
		    continue;
		}
		bool trimstatus = BRNode_isTrimmed(list[i], uv, &v);
		if (v >= 0.0) {
		    if (closest && *closest == NULL) {
			currHeight = v;
			currTrimStatus = trimstatus;
			if (closest) 
			    **closest = list[i];
		    } else if (v < currHeight) {
			currHeight = v;
			currTrimStatus = trimstatus;
			if (closest) 
			    **closest = list[i];
		    }
		} else {
		    dist = fabs(v);
		    if (!underTrim) {
			underTrim = true;
			udist = dist;
			*uclosest = list[i];
		    } else {
		    if (udist > dist) 
		    	udist = dist;
			*uclosest = list[i].;
		    }
		}
	    }

	    if (closest && *closest == NULL) {
			if (verticalTrim) {
			    *closesttrim = vdist;
			    if (closest) 
					*closest = vclosest; 
			}
			if ((underTrim) && (!verticalTrim || (udist < *closesttrim))) {
			    *closesttrim = udist;
			    if (closest) 
					*closest = uclosest;
			}
			return true;
	    } else {
			*closesttrim = currHeight;
			if ((verticalTrim) && (vdist < *closesttrim)) {
			    *closesttrim = vdist;
			    if (closest) {
				*closest = vclosest;
			    }
			}
			if ((underTrim) && (udist < *closesttrim)) {
			    *closesttrim = udist;
			    if (closest) {
				*closest = uclosest;
			    }
			}
			return currTrimStatus;
		    }
	}
    } else {
	if (sbv.m_trimmed)
	    return true;
	else 
	    return false;
    }
}




bool intersectedBy(struct BBNode sbv, double3 r_pt_, double3 r_dir_)
{
    double tnear = -DBL_MAX;
    double tfar = DBL_MAX;
    bool untrimmedresult = true;
    double r_pt[3]={r_pt_.x, r_pt_.y, r_pt_.z};
    double r_dir[3]={r_dir_.x, r_dir_.y, r_dir_.z};
    for (int i = 0; i < 3; i++) {
	if (UNLIKELY(NEAR_ZERO(r_dir[i]))) {
	    if (r_pt[i] < sbv.m_min[i] || r_pt[i] > sbv.m_max[i])
			untrimmedresult = false;
	    
	} else {
	    double t1 = (sbv.m_min[i] - r_pt[i]) / r_dir[i];
	    double t2 = (sbv.m_max[i] - r_pt[i]) / r_dir[i];
	    if (t1 > t2) {
		double tmp = t1;    /* swap */
		t1 = t2;
		t2 = tmp;
	    }

	    tnear = fmax(tnear, t1);
	    tfar = fmin(tfar, t2);

	    if (tnear > tfar) /* box is missed */
			untrimmedresult = false;
	}
    }
    if (m_stl->m_children.empty())
		return !sbv.m_trimmed && untrimmedresult;
    else
		return untrimmedresult;
    
}



int intersectsHierarchy(struct BBNode sbv, double3 r_pt, double3 r_dir, struct BBNode* results_opt)
{
    bool intersects = intersectedBy(sbv, r_pt, r_dir);
    if (intersects && sbv.m_stl->m_children.empty())
		results_opt.pushback(this);
    else if (intersects) {
		for (int i = 0; i < sbv.m_stl->m_children.size(); i++) 
		    intersectsHierarchy(sbv.m_stl->m_children[i], r_pt, r_dir, results_opt);
    }
    
}



void utah_pushBack(struct BBNode* sbv, double3* uv)
{
    double t0, t1;
    int i = sbv->m_u[0] < sbv->m_u[1] ? 0 : 1;

    t0 = sbv->m_u[i];
    t1 = sbv->m_u[1 - i];
    if (*uv.x < t0)
		*uv.x = t0;
    else if (*uv.x > t1) 
		*uv.x = t1;
    
    i = sbv->m_v[0] < sbv->m_v[1] ? 0 : 1;
    t0 = sbv->m_v[i];
    t1 = sbv->m_v[1 - i];
    if (*uv.y < t0) 
		*uv.y = t0;
    else if (*uv.y > t1) 
		*uv.y = t1;
    
}


int
utah_newton_solver(struct BBNode* sbv, struct ON_Surface* surf, double3 r_pt, double3 r_dir, double3* ouv, double* t, double3* N, bool* converged, double3* suv, int count, int iu, int iv)
{
    int i = 0;
    int intersects = 0;
    double j11 = 0.0;
    double j12 = 0.0;
    double j21 = 0.0;
    double j22 = 0.0;
    double f = 0.0;
    double g = 0.0;
    double rootdist = 0.0;
    double oldrootdist = 0.0;
    double J = 0.0;
    double invdetJ = 0.0;
    double du = 0.0;
    double dv = 0.0;
    double cdu = 0.0;
    double cdv = 0.0;

    double3 p1, p2;
    double p1d = 0.0, p2d = 0.0;
    int errantcount = 0;

    if (fabs(r_dir.x) > fabs(r_dir.y) && fabs(r_dir.x) > fabs(r_dir.z))
		p1 = (double3)(r_dir.y, -r_dir.x, 0);
    else
		p1 = (double3)(0, r_dir.z, -r_dir.y);
    p1 = normalize(p1);

    p2 = cross(p1, r_dir);

    p1d = -dot(p1, r_pt);
    p2d = -dot(p2, r_pt);

    double3 S = (double3)(0.0);
    double3 Su = (double3)(0.0);
    double3 Sv = (double3)(0.0);

    double3 uv = (double3)(0.0);
    double3 puv = (double3)(0.0);

    uv.x = suv->x;
    uv.y = suv->y;

    double3 uv0 = uv;
    surface_Ev1Der(uv.x, uv.y, &S, &Su, &Sv, 0, 0);
    
    f = dot(S, p1) + p1d;
    g = dot(S, p2) + p2d;
    rootdist = fabs(f) + fabs(g);

    for (i = 0; i < BREP_MAX_ITERATIONS; i++) {
    j11 = dot(Su, p1);
    j21 = dot(Su, p2);
    j12 = dot(Sv, p1);
    j22 = dot(Sv, p2);

	J = (j11 * j22 - j12 * j21);

	if (NEAR_ZERO(J, BREP_INTERSECTION_ROOT_EPSILON)) {
	    // perform jittered perturbation in parametric domain....
	    uv.x = uv.x + .1 * drand48() * (uv0.x - uv.x);
	    uv.y = uv.y + .1 * drand48() * (uv0.y - uv.y);
	    continue;
	}

	invdetJ = 1.0 / J;

	if ((iu != -1) && (iv != -1)) {
	    du = -invdetJ * (j22 * f - j12 * g);
	    dv = -invdetJ * (j11 * g - j21 * f);

	    if (i == 0) {
		if (((iu == 0) && (du < 0.0)) || ((iu == 1) && (du > 0.0)))
		    return intersects; //head out of U bounds
		if (((iv == 0) && (dv < 0.0)) || ((iv == 1) && (dv > 0.0)))
		    return intersects; //head out of V bounds
	    }
	}

	du = invdetJ * (j22 * f - j12 * g);
	dv = invdetJ * (j11 * g - j21 * f);

	if (i!=0) {
	    int sgnd = (du > 0) - (du < 0);
	    int sgncd = (cdu > 0) - (cdu < 0);
	    if ((sgnd != sgncd) && (fabs(du) > fabs(cdu))) 
			du = sgnd * 0.75 * fabs(cdu);
	    
	    sgnd = (dv > 0) - (dv < 0);
	    sgncd = (cdv > 0) - (cdv < 0);
	    if ((sgnd != sgncd) && (fabs(dv) > fabs(cdv))) 
			dv = sgnd * 0.75 * fabs(cdv);
	}	    
	
	cdu = du;
	cdv = dv;
	
	puv.x = uv.x;
	puv.y = uv.y;

	uv.x -= du;
	uv.y -= dv;

	utah_pushBack(sbv, &uv);

	surface_Ev1Der(uv.x, uv.y, &S, &Su, &Sv, 0, 0);
    f = (S * p1) + p1d;
    g = (S * p2) + p2d;
	oldrootdist = rootdist;
	rootdist = fabs(f) + fabs(g);
	int halve_count = 0;

	/* FIXME: all constants should be documented, why this
	 * value?  what's the sensitivity/impact?
	 */
	while ((halve_count++ < 3) && (oldrootdist < rootdist)) {
	    // divide current UV step
	    uv.x = (puv.x + uv.x) / 2.0;
	    uv.y = (puv.y + uv.y) / 2.0;

	    utah_pushBack(sbv, &uv);

	    surface_Ev1Der(uv.x, uv.y, &S, &Su, &Sv, 0, 0);
	    
	    f = (S * p1) + p1d;
    	g = (S * p2) + p2d;
	    rootdist = fabs(f) + fabs(g);
	}

	if (oldrootdist <= rootdist) {

	    /* FIXME: all constants should be documented. why this
	     * value? must it coincide with the constant in the
	     * preceding loop?
	     */
	    if (errantcount > 3) 
			return intersects;
	    else 
			errantcount++;
	    
	}

	if (rootdist < ROOT_TOL) {
	    int ulow = (sbv->m_u[0] <= sbv->m_u[1]) ? 0 : 1;
	    int vlow = (sbv->m_v[0] <= sbv->m_v[1]) ? 0 : 1;
	    if ((sbv->m_u[ulow] - VUNITIZE_TOL < uv.x && uv.x < sbv->m_u[1 - ulow] + VUNITIZE_TOL) &&
		(sbv->m_v[vlow] - VUNITIZE_TOL < uv.y && uv.y < sbv->m_v[1 - vlow] + VUNITIZE_TOL)) {
		bool new_point = true;
		for (int j = 0; j < count; j++) {
		    if (NEAR_EQUAL(uv.x, ouv[j].x, VUNITIZE_TOL) && NEAR_EQUAL(uv.y, ouv[j].y, VUNITIZE_TOL)) 
				new_point = false;
		    
		}
		if (new_point) {
		    t[count] = dot(r_dir, (S-r_pt))/dot(r_dir, r_dir);
		    N[count] = cross(Su, Sv);
		    N[count] = normalize(N[count]);
		    ouv[count].x = uv.x;
		    ouv[count].y = uv.y;
		    intersects++;
		    *converged = true;
		}
	    }
	    return intersects;
	}
    }
    return intersects;
}


int
utah_newton_4corner_solver(struct BBNode* sbv, struct ON_Surface* surf, double3 r_pt, double3 r_dir, double3* ouv, double* t, double3* N, bool* converged, int docorners)
{
    int intersects = 0;
    *converged = false;
    if (docorners) {
	for (int iu = 0; iu < 2; iu++) {
	    for (int iv = 0; iv < 2; iv++) {
		double3 uv = (double3)(0.0);
		uv.x = sbv->m_u[iu];
		uv.y = sbv->m_v[iv];
		intersects += utah_newton_solver(sbv, surf, r_pt, r_dir, ouv, t, N, converged, &uv, intersects, iu, iv);
	    }
	}
    }

    double3 uv = (double3)(0.0);
    uv.x = 0.5 * (sbv->m_u[0] + sbv->m_u[1]);
    uv.y = 0.5 * (sbv->m_v[0] + sbv->m_v[1]);
    intersects += utah_newton_solver(sbv, surf, r_pt, r_dir, ouv, t, N, converged, &uv, intersects, -1, -1);
    return intersects;
}


int
utah_brep_intersect(struct BBNode* sbv, struct ON_BrepFace* face, struct ON_Surface* surf, double3* uv, double3 r_pt, double3 r_dir, struct brep_hit *hits, int s)
{
	struct brep_hit *list;
	list = &hits[0];
    double3 N[MAX_BREP_SUBDIVISION_INTERSECTS];
    double t[MAX_BREP_SUBDIVISION_INTERSECTS];
    double3 ouv[MAX_BREP_SUBDIVISION_INTERSECTS];
    int found = BREP_INTERSECT_ROOT_DIVERGED;
    bool converged = false;
    int numhits;

    double grazing_float = dot(sbv->m_normal, r_dir);

    if (fabs(grazing_float) < 0.2) 
		numhits = utah_newton_4corner_solver(sbv, surf, r_pt, r_dir, ouv, t, N, &converged, 1);
    else 
		numhits = utah_newton_4corner_solver(sbv, surf, r_pt, r_dir, ouv, t, N, &converged, 0);
    
    if (converged) {
	for (int i = 0; i < numhits; i++) {
	    double closesttrim;
	    struct BRNode* trimBR = NULL;
	    int trim_status = BBNode_isTrimmed(*sbv, ouv[i], &trimBR, &closesttrim, BREP_EDGE_MISS_TOLERANCE);
	    if (trim_status != 1) {
			double3 _pt;
			double3 _norm = vload3(0, N[i]);
			double3 vpt;
			double3 vnorm;
			_pt = r_pt + (r_dir * t[i]);
			vpt = vload3(0, _pt);
			if (face->m_bRev) 
				_norm = -vload3(0, _norm);
			
			vnorm = _norm;
			uv.x = ouv[i].x;
			uv.y = ouv[i].y;

			list[s].face = face;
			list[s].dist = t[i];
			list[s].oob = false;
			list[s].origin = vload3(0, r_pt);
			list[s].point = vload3(0, vpt);
			list[s].normal = vload3(0, vnorm);
			list[s].sbv = sbv;
			list[s].uv = vload3(0, uv);

			list[s].trimmed = false;
			if (trimBR != NULL) 
			    list[s].m_adj_face_index = trimBR->m_adj_face_index;
			else 
			    list[s].m_adj_face_index = -99;
			
			if (fabs(closesttrim) < BREP_EDGE_MISS_TOLERANCE) {
			    list[s].closeToEdge = true;
			    list[s].hit = NEAR_HIT;
			} else {
			    list[s].closeToEdge = false;
			    list[s].hit = CLEAN_HIT;
			}
			if (dot(r_dir, vnorm) < 0.0)
			    list[s].direction = ENTERING;
			else
			    list[s].direction = LEAVING;
			
			found = BREP_INTERSECT_FOUND;
	    } else if (fabs(closesttrim) < BREP_EDGE_MISS_TOLERANCE) {
			double3 _pt;
			double3 _norm = vload3(0, N[i]);
			double3 vpt;
			double3 vnorm;
			_pt = r_pt + (r_dir * t[i]);
			vpt = vload3(0, _pt);
			if (face->m_bRev) 
				_norm = -vload3(0, _norm);
			vnorm = _norm;
			uv.x = ouv[i].x;
			uv.y = ouv[i].y;

			list[s].face = face;
			list[s].dist = t[i];
			list[s].oob = false;
			list[s].origin = vload3(0, r_pt);
			list[s].point = vload3(0, vpt);
			list[s].normal = vload3(0, vnorm);
			list[s].uv = vload3(0, uv);
			list[s].trimmed = true;
			list[s].closeToEdge = true;
			list[s].sbv = sbv;
			list[s].hit = NEAR_MISS;

			if (trimBR != NULL) 
			    list[s].m_adj_face_index = trimBR->m_adj_face_index;
			else 
			    list[s].m_adj_face_index = -99;
			
			if (dot(r_dir, vnorm) < 0.0)
			    list[s].direction = ENTERING;
			else
			    list[s].direction = LEAVING;
			
			found = BREP_INTERSECT_FOUND;
	    }
	}
    }
    return found;
}






int containsNearMiss(struct brep_hit *hits, int total_hits)
{
	struct brep_hit *list;
    list=&hits[0];
    for (int i=0; i<total_hits;i++) {
	if (list[i].hit == NEAR_MISS && list[i].useful) 
	    return 1;
    }
    return 0;
}


int containsNearHit(struct brep_hit *hits, int total_hits)
{
	struct brep_hit *list;
    list=&hits[0];
    for (int i=0; i<total_hits;i++) {
	if (list[i].hit == NEAR_HIT && list[i].useful) 
	    return 1;
    }
    return 0;
}

void deleteNode(struct brep_hit *hits, int curr, int *front, int *back, int *true_hits)
{
	struct brep_hit *list;
    list=&hits[0];
    list[curr].useful=0;
    (*true_hits)--;
    if(curr==(*front))
    {
    	while(list[(*front)].useful==0)
    		(*front)++;
    }
    
    if(curr==(*back))
    {
    	while(list[(*back)].useful==0)
    		(*back)--;
    }
}

int findPrev(struct brep_hit *hits, int ctr)
{
	struct brep_hit *list;
    list=&hits[0];
    while(list[ctr].useful==0)
    	ctr--;
    return ctr;
}

int findNext(struct brep_hit *hits, int ctr)
{
	struct brep_hit *list;
    list=&hits[0];
    while(list[ctr].useful==0)
    	ctr++;
    return ctr;
}

/**
 * Intersect a ray with a brep.  If an intersection occurs, a struct
 * seg will be acquired and filled in.
 *
 * Returns -
 * 0 MISS
 * >0 HIT
 */



int brep_shot(RESULT_TYPE *res, const double3 r_pt, const double3 r_dir, const uint idx, global const struct brep_specific *bs)
{
    
    int total_hits = 0; int true_hits;
    int curr;
    int front, back;
    int prev, next;
    struct brep_hit *list;
    struct brep_hit *all_hits; // record all hits
    int len;

    if (!bs)
		return 0;

    /* First, test for intersections between the Surface Tree
     * hierarchy and the ray - if one or more leaf nodes are
     * intersected, there is potentially a hit and more evaluation is
     * needed.  Otherwise, return a miss.
     */
    struct BBNode* inters;
    len = intersectsHierarchy(bs->bvh, r_pt, r_dir, inters);
    if (len==0) return 0; // MISS

    for (int i=0;i <len;i++) {
	     struct BBNode* sbv = &inters[i];
	     struct ON_BrepFace* f = &sbv->m_face;
	     struct ON_Surface* surf = f->SurfaceOf();
	     double3 uv = (0.5 * (sbv->m_u[0] + sbv->m_u[1]), 0.5 * (sbv->m_v[0] + sbv->m_v[1]), 0.0);
	     if(utah_brep_intersect(sbv, f, surf, uv, r_pt, r_dir, all_hits, total_hits) == BREP_INTERSECT_FOUND)
		      total_hits++;
    }

    /* sort the hits */
    primitive_hitsort(all_hits, total_hits);
  
    list=&all_hits[0];
    true_hits = total_hits;
    front=0;
    back=total_hits-1;

    if ((true_hits > 1) && containsNearMiss(all_hits, total_hits)) { 
    	curr=front;
    	while (curr <= back) {

	    if (list[curr].hit == NEAR_MISS && list[curr].useful) {
		    if (curr != front) {  
			     prev = findPrev(all_hits, curr-1);
		       if ((list[prev].hit != NEAR_MISS) && (list[prev].direction == list[curr].direction)) {
			       //remove current miss
		         deleteNode(all_hits, curr, &front, &back, &true_hits);
			       curr=front; //rewind and start again
			       continue;
		    }
		}
		
		if (curr != back) { 
			next = findNext(all_hits, curr+1);
		    if ((list[next].hit != NEAR_MISS) && (list[next].direction == list[curr].direction)) {
			//remove current miss
			deleteNode(all_hits, curr, &front, &back, &true_hits);
			curr=front; //rewind and start again   
			continue;
		    }
		}
	    }
	    curr++;
	}

	curr=front;
	// check for crack hits between adjacent faces
	while (curr <= back) {
	    if (curr != front) {
		if ((list[curr].hit == NEAR_MISS) && list[curr].useful) {
			prev = findPrev(all_hits, curr-1);
		    if (list[prev].hit == NEAR_MISS) { // two near misses in a row
			if (list[prev].m_adj_face_index == list[curr].face.m_face_index) {
			    if (list[prev].direction == list[curr].direction) {
				//remove current miss
				list[prev].hit = CRACK_HIT;
				deleteNode(all_hits, curr, &front, &back, &true_hits);
				curr++;
				continue;
			    } else {
				//remove both edge near misses
				deleteNode(all_hits, prev, &front, &back, &true_hits);
				deleteNode(all_hits, curr, &front, &back, &true_hits);
				curr++;
				continue;
			    }
			} else {
			    // not adjacent faces so remove first miss
				deleteNode(all_hits, prev, &front, &back, &true_hits);
			}
		    }
		}
	    }
	    curr++;
	}
	// check for CH double enter or double leave between adjacent faces(represents overlapping faces)

	curr=front;
	while (curr <= back) {

	    if (list[curr].hit == CLEAN_HIT && list[curr].useful) {
		if (curr != front) {
			prev = findPrev(all_hits, curr-1);
		    if ((list[prev].hit == CLEAN_HIT) && (list[prev].direction == list[curr].direction) && (list[prev].face.m_face_index == list[curr].m_adj_face_index)) {
			// if "entering" remove first hit if "exiting" remove second hit
			// until we get good solids with known normal directions assume
			// first hit direction is "entering" todo check solid status and normals

			if (list[front].direction == list[curr].direction)  // assume "entering"
			    {deleteNode(all_hits, prev, &front, &back, &true_hits);}
			else  // assume "exiting"
			    {deleteNode(all_hits, curr, &front, &back, &true_hits); curr++;}
			
			continue;
		    }
		}
	    }
	    curr++;
	}

	if ((true_hits>0) && ((true_hits % 2) != 0) && (list[back].hit == NEAR_MISS)) 
			deleteNode(all_hits, back, &front, &back, &true_hits);
			
	if ((true_hits>0) && ((true_hits % 2) != 0) && (list[front].hit == NEAR_MISS))
		    deleteNode(all_hits, front, &front, &back, &true_hits);

    }



    ///////////// handle near hit
    if ((true_hits > 1) && containsNearHit(all_hits, total_hits)) { //&& ((true_hits % 2) != 0)) {

	curr=front;
	while (curr <= back) {
	    if (list[curr].hit == NEAR_HIT && list[curr].useful) {
		if (curr != front) {
			prev = findPrev(all_hits, curr-1);
		    if ((list[prev].hit != NEAR_HIT) && (list[prev].direction == list[curr].direction)) {
			//remove current miss
		    deleteNode(all_hits, curr, &front, &back, &true_hits);
		    }
		}

		else if (curr != back) {
			next = findNext(all_hits, curr+1);
		    if ((list[next].hit != NEAR_HIT) && (list[next].direction == list[curr].direction)) {
			//remove current miss
		    deleteNode(all_hits, curr, &front, &back, &true_hits);
		    }
		}
	    }
	    curr++;
	}

	curr=front;
	while (curr <= back) {
	    if (list[curr].hit == NEAR_HIT && list[curr].useful) {
		if (curr != front) {
			prev = findPrev(all_hits, curr-1);
		    if ((list[prev].hit == NEAR_HIT) && (list[prev].direction == list[curr].direction)) {
			//remove current near hit
			list[prev].hit = CRACK_HIT;
		    deleteNode(all_hits, curr, &front, &back, &true_hits);
		    }
		}
	    }
	    curr++;
	}
    }


    if (true_hits>0) {
	// remove grazing hits with with normal to ray dot less than BREP_GRAZING_DOT_TOL (>= 89.999 degrees obliq)
	for (curr=front; curr<=back; curr++) {
	    if(list[curr].useful)
	    {
		    if ((list[curr].trimmed && !list[curr].closeToEdge) || list[curr].oob || NEAR_ZERO(dot(list[curr].normal, r_dir), BREP_GRAZING_DOT_TOL)) {
				// remove what we were removing earlier
				deleteNode(all_hits, curr, &front, &back, &true_hits);
				curr--;
		    }
		}
	}
    }


    if (true_hits>0) {
	// we should have "valid" points now, remove duplicates or grazes(same point with in/out sign change)
	prev=0;
	curr=1;
	while (curr <= back) {
	    if(list[curr].useful)
	    {	
	    if ((NEAR_ZERO(list[curr].dist - list[prev].dist, BREP_SAME_double3OLERANCE)) && list[prev].useful) {
			double prevDot = dot(list[prev].normal, r_dir);
			double iDot = dot(list[curr].normal, r_dir);

			if ((prevDot*iDot)<0) {
			    // delete them both
		  		deleteNode(all_hits, prev, &front, &back, &true_hits);
		   		deleteNode(all_hits, curr, &front, &back, &true_hits);
			    curr++;
			    prev = curr;
			} else {
			    // just delete the second
		  		deleteNode(all_hits, curr, &front, &back, &true_hits);
			}
	    } else 
			prev = curr;
		}
		curr++;
	}
    }

    // remove multiple "INs" in a row assume prev "IN" is the actual entering hit, for
    // multiple "OUTs" in a row assume first "OUT" is the actual exiting hit, remove unused
    // "INs/OUTs" from hit list.
    //if (true_hits>0 && ((true_hits % 2) != 0)) {
    if (true_hits>0) {
	// we should have "valid" points now, remove duplicates or grazes
	int prev=0;
	curr=1;
	int entering = 1;

	while (curr <= back) {
	    double prevDot = dot(list[prev].normal, r_dir);
	    double iDot = dot(list[curr].normal, r_dir);		    

		if(list[curr].useful)
		{
		    if (curr == front) 
				entering = iDot;
		    
		    if ((prevDot*iDot >0) && list[prev].useful) {
				if (iDot*entering >0) {
				    deleteNode(all_hits, prev, &front, &back, &true_hits);
				    prev = curr;
				} else  //exiting
				    deleteNode(all_hits, curr, &front, &back, &true_hits);
			    
		    } else 
				prev = curr;
		}
	curr++;
	}
    }


    if ((true_hits > 1) && ((true_hits % 2) != 0)) {
		double firstDot = dot(list[front].normal, r_dir);
		double prevDot = dot(list[back].normal, r_dir);
		if ((firstDot*prevDot)>0) 
			    deleteNode(all_hits, back, &front, &back, &true_hits);
	
    }

    if ((true_hits > 1) && (true_hits % 2 == 0)){

    	struct hit hits[2];
	    // take each pair as a segment
	    for (curr=front; curr<=back;curr++) {
	    	
			if(list[curr].useful==0)
				continue;

			struct brep_hit in = list[curr];
			curr++;

			while(list[curr].useful==0)
				curr++; 
				
			struct brep_hit out = list[curr];
			
			hits[0].hit_point = vload3(0, in.point);
			hits[0].hit_normal = vload3(0, in.normal);
			hits[0].hit_dist = in.dist;
			hits[0].hit_surfno = in.face.m_face_index;
			hits[0].hit_vpriv = vload3(0, in.uv);

			hits[1].hit_point = vload3(0, out.point);
			hits[1].hit_normal = vload3(0, out.normal);
			hits[1].hit_dist = out.dist;
			hits[1].hit_surfno = out.face.m_face_index;
			hits[1].hit_vpriv = vload3(0, out.uv);

			do_segp(res, idx, &hits[0], &hits[1]);
	    }
	    return true_hits;
	} 

    return 0; // MISS
}