// TO BE INCLUDED IN rt.cl

/**
 * check if a ray passes within a bounding sphere
 */
int rt_in_sph(const double3 r_pt, const double3 r_dir, double3 center, double radius_sq)
{
    double3 toCenter;
    double3 toPCA;
    double dist_sq;

    toCenter = center - r_pt;
    toPCA = cross(toCenter, r_dir);
    dist_sq = dot(toPCA, toPCA);

    if (dist_sq <= radius_sq) 
	   return 1;
    else 
	   return 0;
    
}
