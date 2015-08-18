#ifndef mode_h__
#define mode_h__

extern void findmode(
    //mode-specific variables
    long order, double modebounces[],
    //initial position on PSS
    double theta0, double sinchi0,
    //root-finding algorithm parameters
    double rtol, long maxiter,
    //Cavity parameters
    void *bnd, double (*rfunc_p)(void *bnd,double theta),
    void (*rsys_p)(void *bnd, double theta,double results[]),
    void *idx, double (*nfunc_p)(void *idx, double r, double theta),
    void (*nderiv_p)(void *idx, double r, double theta,double results[]));

#endif