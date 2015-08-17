#ifndef ray_h__
#define ray_h__

extern void rayevolve(
    //Storage arrays
    double raypath_r[], double raypath_theta[], int bounces_indices[],
    double bounces_chi[], int lengths[],
    //Initial conditions
    double r0, double theta0, double pr0, double ptheta0,
    //Simulation parameters
    double tmax, int bouncemax, double reltol, double abstol,
    //Cavity parameters
    void *bnd, double (*rfunc_p)(void *bnd,double theta),
    void (*rsys_p)(void *bnd, double theta,double results[]),
    void *idx, double (*nfunc_p)(void *idx, double r, double theta),
    void (*nderiv_p)(void *idx, double r, double theta,double results[]));

#endif