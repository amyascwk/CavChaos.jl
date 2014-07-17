#ifndef ray_h__
#define ray_h__

extern void rayevolve(
    double raypath0[], double raypath1[], int bounces0[], double bounces1[],
    int lengths[], double r0, double theta0, double pr0, double ptheta0,
    double tmax, int bouncemax, double reltol, double abstol,
    void *bndthunk,
    double (*rfunc)(void *bndthunk,double theta),
    void (*rsys)(void *bndthunk, double theta,double results[]),
    void *idxthunk,
    double (*nfunc)(void *idxthunk, double r, double theta),
    void (*nsys)(void *idxthunk, double r, double theta,double results[]));

#endif