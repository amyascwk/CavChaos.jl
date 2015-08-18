#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

#define pi M_PI


//#############################################################################
//#############################################################################
//Ray path solver

//Data structure for pointers to julia cavity data and functions
typedef struct {
    //pointer to Boundary object
    void *bnd;
    //pointer to radius function
    double (*rfunc_p)(void *bnd,double theta);
    //pointer to (mutating) radius and normal vector angle function
    void (*rsys_p)(void *bnd, double theta, double results[]);
    //pointer to RefractiveIndex object
    void *idx;
    //pointer to refractive index value function
    double (*nfunc_p)(void *idx, double r, double theta);
    //pointer to (mutating) refractive index value and derivative function
    void (*nderiv_p)(void *idx, double r, double theta, double results[]);
} cavinfo;


//#############################################################################
//Functions for ODE

//ODE derivatives function
int odefunc(double t, const double y[], double f[], void *params){
    //Get Julia pointers
    cavinfo *cip = (cavinfo *)params;
    
    //Calculate n and its derivatives and store in results array
    double results[3];
    (*(*cip).nderiv_p)((*cip).idx,y[0],y[1],results);
    //results[0] = n, results[1] = dn/dr, results[2] = dn/dtheta
    
    //Debugging purposes
    //printf("n = %.5f, dr_n = %.5f, dtheta_n = %.5f\n",n,dr_n,dtheta_n);
    
    //Calculate and store derivative of ODE coordinate vector
    //Equations are:
    //r' = p_r/n
    //theta' = p_theta/(r^2*n)
    //p_r' = p_theta^2/(r^3*n) + dn/dr
    //p_theta' = dn/dtheta
    f[0] = y[2]/results[0];
    f[1] = y[3]/(y[0]*y[0]*results[0]);
    f[2] = y[3]*y[3]/(y[0]*y[0]*y[0]*results[0]) + results[1];
    f[3] = results[2];
    
    return GSL_SUCCESS;
}


//ODE Jacobian function
int odejac(double t, const double y[], double *dfdy, double dfdt[], void *params){
    //RK8PD does not require the Jacobian
    return GSL_SUCCESS;
}


//#############################################################################
//Solver functions

//Ray reflection
//Re-compute a ray coordinate vector to simulate a bounce, using the ODE coordinate vectors S0 and S immediately before and after the bounce, where S0 = (r0,theta0,pr0,ptheta0), and S = (r,theta,pr,ptheta).
double raybounce(cavinfo *cip, double S0[], double S[]){
    
    //Binary search for intersection of trajectory with cavity boundary
    //Assume that light trajectory between (r0,theta0) and (r,theta) is a straight line
    //(fastest way to interpolate, after all this is already within 1 stepsize of time)
    const double x0 = S0[0]*cos(S0[1]), y0 = S0[0]*sin(S0[1]);
    const double x1 = S[0]*cos(S[1]), y1 = S[0]*sin(S[1]);
    double uA = 0.0, uB = 1.0; //bounds of the binary search
    double xC,yC,rC,thetaC,RC,uC;
    do{
        uC = 0.5*(uA+uB); //get middle point
        xC = (1-uC)*x0 + uC*x1; yC = (1-uC)*y0 + uC*y1;
        rC = hypot(xC,yC); thetaC = atan2(yC,xC);
        RC = (*(*cip).rfunc_p)((*cip).bnd,thetaC);
        //Change boundary
        if(rC > RC) uB = uC;
        else uA = uC;
    } while(fabs(rC-RC) > 1e-12);
    
    //Get angle of incidence and store in results array
    double results[2];
    (*(*cip).rsys_p)((*cip).bnd,thetaC,results);
    //results[0] = rC, results[1] = alpha
    //Use linear interpolation to find ray angle phi at intersection
    const double phi0 = S0[1] + atan2(S0[3],S0[0]*S0[2]);
    const double phi1 = S[1] + atan2(S[3],S[0]*S[2]);
    const double phiC = (1-uC)*phi0 + uC*phi1;
    //Angle of incidence
    const double chi = phiC - results[1];
    
    //Store interpolated and reflected ODE coordinate vector
    const double phi = pi - chi + results[1]; //ray angle after reflection
    const double n = (*(*cip).nfunc_p)((*cip).idx,rC,thetaC);
    const double pr = n*cos(phi-thetaC);
    const double ptheta = n*rC*sin(phi-thetaC);
    S[0] = rC; S[1] = thetaC; S[2] = pr; S[3] = ptheta;
    
    //Report bounce information
    return chi;
}


//Ray evolution
void rayevolve(
    //Storage arrays
    double raypath_r[], double raypath_theta[], long bounceindices[],
    double bouncepts_chi[], long lengths[],
    //Initial conditions
    double r0, double theta0, double pr0, double ptheta0,
    //Simulation parameters
    double tmax, long bouncemax, double reltol, double abstol,
    //Cavity parameters
    void *bnd, double (*rfunc_p)(void *bnd,double theta),
    void (*rsys_p)(void *bnd, double theta,double results[]),
    void *idx, double (*nfunc_p)(void *idx, double r, double theta),
    void (*nderiv_p)(void *idx, double r, double theta,double results[])){
        
        //Condense input cavity data into cavinfo struct
        cavinfo ci = {bnd,rfunc_p,rsys_p,idx,nfunc_p,nderiv_p};
        
        //Initialize solver
        gsl_odeiv2_system sys = {odefunc,odejac,4,&ci};
        const gsl_odeiv2_step_type *steptype = gsl_odeiv2_step_rk8pd;
        gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(steptype,4);
        gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(abstol,reltol);
        gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(4);
        
        //Initialize parameters
        double t = 0.0, dt = 0.0001;
        double y0[4], y[4] = {r0,theta0,pr0,ptheta0};
        
        //Prepare results record
        long bouncenum = 0, stepnum = 1; //indicates postion to record next
        const long prealloc = 250*ceil(tmax); //length of preallocated raypath array
        raypath_r[0] = y[0]; raypath_theta[0] = y[1];
        
        //Solver loop
        while(t < tmax && bouncenum < bouncemax && stepnum < prealloc){
            //Record initial position
            y0[0] = y[0], y0[1] = y[1], y0[2] = y[2], y0[3] = y[3];
            
            //Run Solver
            int status = 
                gsl_odeiv2_evolve_apply(evolve,control,step,&sys,&t,tmax,&dt,y);
            if(status != GSL_SUCCESS) break;
            
            //Check Hamiltonian once in a while for sanity check
            if(stepnum%1000 == 0){
                double H = y[2]*y[2]+y[3]*y[3]/(y[0]*y[0]) - gsl_pow_2((*ci.nfunc_p)(ci.idx,y[0],y[1]));
                if(H > 1e-9) printf("Warning: Error in Hamiltonian is %.f\n",H);
            }
            
            //Check difference in ray and boundary radial positions
            double dr = y[0] - (*ci.rfunc_p)(ci.bnd,y[1]);
            if(dr > 0){
                //Boundary crossing!
                //get chi and corrected y
                bouncepts_chi[bouncenum] = raybounce(&ci,y0,y);
                //store (Julia's 1-based) index for thetaC values recorded in 
                //raypath_theta array
                bounceindices[bouncenum] = stepnum+1;
                bouncenum += 1;
            }
            
            //Record position 
            raypath_r[stepnum] = y[0]; raypath_theta[stepnum] = y[1];
            stepnum += 1;
            
            //Display progress for debugging
            //printf("t = %.5f, y = [%.5f,%.5f,%.5f,%.5f]\n",t,y[0],y[1],y[2],y[3]);
            
        }
        //Store lengths of arrays
        lengths[0] = stepnum; lengths[1] = bouncenum;
        
        //Free memory
        gsl_odeiv2_step_free(step);
        gsl_odeiv2_control_free(control);
        gsl_odeiv2_evolve_free(evolve);
        
}

