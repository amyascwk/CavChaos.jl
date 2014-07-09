#ifndef ray_h__
#define ray_h__

extern void rayevolve(
    double raypath0[], double raypath1[], double bounces0[], double bounces1[],
    double r0, double theta0, double pr0, double ptheta0,
    double tmax, double reltol, double abstol,
    void *rthunk,
    double (*rfunc)(void *rthunk,double theta),
    void *rsysthunk,
    void (*rsys)(void *rsysthunk, double theta,double results[]),
    void *nthunk,
    double (*nfunc)(void *nthunk, double r, double theta), 
    void *nsysthunk,
    void (*nsys)(void *nsysthunk, double r, double theta,double results[]));

#endif