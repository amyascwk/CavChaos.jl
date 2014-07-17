#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>


//=============================================================================
//Definitions

const double pi = M_PI;

//Data structure for julia callbacks
struct juliafuncs_t {
    void *bndthunk;
    double (*rfunc)(void *bndthunk,double theta);
    void (*rsys)(void *bndthunk, double theta, double results[]);
    void *idxthunk;
    double (*nfunc)(void *idxthunk, double r, double theta);
    void (*nsys)(void *idxthunk, double r, double theta, double results[]);
} jf;

//=============================================================================
//Functions for ODE

//ODE derivatives function
int func(double t, const double y[], double f[], void *params){
    double results[3];
    (*jf.nsys)(jf.idxthunk,y[0],y[1],results);
    
    //printf("n = %.5f, dr_n = %.5f, dtheta_n = %.5f\n",n,dr_n,dtheta_n);
    f[0] = y[2]/results[0]; //r' = p_r/n
    f[1] = y[3]/(y[0]*y[0]*results[0]); //theta' = p_theta/(r^2*n)
    f[2] = y[3]*y[3]/(y[0]*y[0]*y[0]*results[0]) + results[1]; //p_r' = p_theta^2/(r^3*n) + dn/dr
    f[3] = results[2]; //p_theta' = dn/dtheta
    return GSL_SUCCESS;
}

//ODE Jacobian function
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params){
    //RK8PD does not require the Jacobian
    return GSL_SUCCESS;
}

//=============================================================================
//Solver functions

//Ray reflection
double raybounce(double r0,double theta0, double S[]){
    //Re-compute a bounced ray coordinate vector
    // r0, theta0:                     Pre-collision position
    // S[] = {r,theta,p_r,p_theta}:    Post-collision ODE coordinate vector

    //Binary search for intersection of trajectory with cavity boundary
    double x0 = r0*cos(theta0), y0 = r0*sin(theta0);
    double xA = x0, yA = y0;
    double xB = S[0]*cos(S[1]), yB = S[0]*sin(S[1]);
    double xC,yC,rC,thetaC,RC;
    do{
        xC = 0.5*(xA+xB); yC = 0.5*(yA+yB);
        rC = hypot(xC,yC); thetaC = atan2(yC,xC);
        RC = (*jf.rfunc)(jf.bndthunk,thetaC);
        if(rC > RC){
            xB = xC; yB = yC;
        } else {
            xA = xC; yA = yC;
        }
    } while(fabs(rC-RC) > 1e-12);

    //Compute reflected ray
    double results[2];
    (*jf.rsys)(jf.bndthunk,thetaC,results);
    double chi = fmod(atan2(yC-y0,xC-x0)-results[1],2*pi); //angle of incidence
    double phi = pi - chi + results[1]; //angle of reflected light to horizontal
    double n = (*jf.nfunc)(jf.idxthunk,rC,thetaC);
    double pr = n*cos(phi-thetaC);
    double ptheta = n*rC*sin(phi-thetaC);
    
    //Mutate ODE coordinate vector
    S[0] = rC; S[1] = thetaC; S[2] = pr; S[3] = ptheta;
    
    //Report bounce information
    return chi;
}


//Ray traversal
void rayevolve(
    double raypath0[], double raypath1[], int bounces0[], double bounces1[],
    int lengths[], double r0, double theta0, double pr0, double ptheta0,
    double tmax, int bouncemax, double reltol, double abstol,
    void *bndthunk,
    double (*rfunc)(void *bndthunk,double theta),
    void (*rsys)(void *bndthunk, double theta,double results[]),
    void *idxthunk,
    double (*nfunc)(void *idxthunk, double r, double theta),
    void (*nsys)(void *idxthunk, double r, double theta,double results[])){
        
        //Condense input functions into juliafuncs_t struct
        jf.bndthunk = bndthunk;
        jf.rfunc = rfunc;
        jf.rsys = rsys;
        jf.idxthunk = idxthunk;
        jf.nfunc = nfunc;
        jf.nsys = nsys;
        
        //Initialize solver
        int *null = NULL;
        gsl_odeiv2_system sys = {func,jac,4,null};
        const gsl_odeiv2_step_type *steptype = gsl_odeiv2_step_rk8pd;
        gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(steptype,4);
        gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-12,1e-12);
        gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(4);
        
        double t = 0.0, dt = 0.0001;
        double y[4] = {r0,theta0,pr0,ptheta0};
        
        //Prepare results record
        int bouncenum = 0, stepnum = 1; //indicates postion to record next
        const int prealloc = 250*ceil(tmax);
        raypath0[0] = y[0]; raypath1[0] = y[1];
        
        //Solver loop
        while(t < tmax && bouncenum < bouncemax && stepnum < prealloc){
            //Remember initial position
            double r0 = y[0], theta0 = y[1];
            
            //Run Solver
            int status = gsl_odeiv2_evolve_apply(evolve,control,step,&sys,&t,tmax,&dt,y);
            if(status != GSL_SUCCESS) break;
            
            //Check Hamiltonian
            if(stepnum%1000 == 0){
                double H = y[2]*y[2]+y[3]*y[3]/(y[0]*y[0]) - gsl_pow_2((*jf.nfunc)(jf.idxthunk,y[0],y[1]));
                if(H > 1e-9) printf("Warning: Error in Hamiltonian is %.f\n",H);
            }
            
            //Check for boundary crossing
            double dr = y[0] - (*jf.rfunc)(jf.bndthunk,y[1]);
            if(dr > 0){
                //get chi and alter y
                bounces1[bouncenum] = raybounce(r0,theta0,y); 
                //store (Julia's 1-based) index for thetaC, thetaC values can be 
                //obtained from raypath.
                bounces0[bouncenum] = stepnum+1;
                bouncenum += 1;
            }
            
            //Record position 
            raypath0[stepnum] = y[0]; raypath1[stepnum] = y[1];
            stepnum += 1;
            
            //Display progress
            //printf("t = %.5f, y = [%.5f,%.5f,%.5f,%.5f]\n",t,y[0],y[1],y[2],y[3]);
            
        }
        //Store lengths of arrays
        lengths[0] = stepnum; lengths[1] = bouncenum;
        
        //Free memory
        gsl_odeiv2_step_free(step);
        gsl_odeiv2_control_free(control);
        gsl_odeiv2_evolve_free(evolve);
        
}
