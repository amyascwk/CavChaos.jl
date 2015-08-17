#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "ray.h"

#define pi M_PI


//#############################################################################
//#############################################################################
//Mode finder

//Data structure for julia cavity data and mode parameters
typedef struct {
    //Julia cavity info and functions -----------------------------------------
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
    
    //mode-specific parameters ------------------------------------------------
    int order; //number of bounces in a period
    double rtol; //stopping criterion tolerance
    
    //preallocated arrays -----------------------------------------------------
    double *raypath_r;
    double *raypath_theta;
    int *bounces_indices;
    double *bounces_chi;
    int *lengths;
    double *modebounces;
} modeinfo;


//#############################################################################
//Function for rootfinding
int rootfunc(const gsl_vector *x, void *params, gsl_vector *f){
    
    //params is pointer to modeinfo struct containing both cavity info and 
    //simulation parameters
    modeinfo *mip = (modeinfo *)params;
    const int order = (*mip).order;
    double *raypath_r = (double *)(*mip).raypath_r;
    double *raypath_theta = (double *)(*mip).raypath_theta;
    int *bounces_indices = (int *)(*mip).bounces_indices;
    double *bounces_chi = (double *)(*mip).bounces_chi;
    int *lengths = (int *)(*mip).lengths;
    
    //Extract initial condition data (print for debugging)
    const double theta0 = gsl_vector_get(x,0), chi0 = gsl_vector_get(x,1);
    //printf("theta0 = %.8f, chi0 = %.8f\n",theta0,chi0);
    
    //Check for sensible values
    if(fabs(chi0) > pi/2){
        //Out of domain
        return GSL_EDOM;
    }
    
    //Setup initial ODE coordinate vector components
    double results[2];
    (*(*mip).rsys_p)((*mip).bnd,theta0,results);
    //result[0] = r0, results[1] = normang0
    const double n = (*(*mip).nfunc_p)((*mip).idx,results[0],theta0);
    const double phi0 = results[1] + chi0;
    const double pr0 = n*cos(phi0-theta0);
    const double ptheta0 = n*results[0]*sin(phi0-theta0);
    
    //Compute ray
    rayevolve(
        //Storage arrays
        raypath_r,raypath_theta,bounces_indices,bounces_chi,lengths,
        //Initial conditions
        results[0],theta0,pr0,ptheta0,
        //Simulation parameters
        200.0,order+1,1e-12,1e-12,
        //Cavity parameters
        (*mip).bnd,(*mip).rfunc_p,(*mip).rsys_p,
        (*mip).idx,(*mip).nfunc_p,(*mip).nderiv_p);
    
    //Store difference in initial and final bounce locations (print for debugging)
    const double dtheta = fmod(raypath_theta[bounces_indices[order]-1] - raypath_theta[bounces_indices[0]-1] + pi, 2*pi) - pi;
    //Equivalent to sin(bounces_chi[order]) - sin(bounces_chi[0]);
    //but more accurate for small differences
    const double dsinchi = 2*sin(0.5*(bounces_chi[order]-bounces_chi[0]))*
        cos(0.5*(bounces_chi[order]+bounces_chi[0]));
    gsl_vector_set(f,0,dtheta/pi);
    gsl_vector_set(f,1,dsinchi);
    //printf("dtheta = %.8f, dsinchi = %.8f\n",dtheta,dsinchi);
    
    //Test stopping criterion
    if(gsl_multiroot_test_residual(f,(*mip).rtol) != GSL_CONTINUE){
        //Record mode bounce data in results array (to avoid running simulation again)
        double *modebounces = (double *)(*mip).modebounces;
        //theta values stored in first half of array (first column after reshaping)
        //sinchi values stored in second half of array (second column after reshaping)
        int i;
        for(i=0; i<order; i+=1){
            modebounces[i] = raypath_theta[bounces_indices[i]-1];
            modebounces[order+i] = sin(bounces_chi[i]);
        }
    }
    
    return GSL_SUCCESS;
}


//#############################################################################
//Solver function
int findmode(
    //mode-specific variables
    int order, double modebounces[],
    //initial position on PSS
    double theta0, double sinchi0,
    //root-finding algorithm parameters
    double rtol, int maxiter,
    //Cavity parameters
    void *bnd, double (*rfunc_p)(void *bnd,double theta),
    void (*rsys_p)(void *bnd, double theta,double results[]),
    void *idx, double (*nfunc_p)(void *idx, double r, double theta),
    void (*nderiv_p)(void *idx, double r, double theta,double results[])){
        
        //Preallocate arrays for reuse
        int bounces_indices[30], lengths[2];
        double raypath_r[50000],raypath_theta[50000],bounces_chi[30];
        
        //Condense input mode data into modeinfo structs
        modeinfo mi = {bnd,rfunc_p,rsys_p,idx,nfunc_p,nderiv_p,order,rtol,
            raypath_r,raypath_theta,bounces_indices,bounces_chi,lengths,modebounces};
        
        //Store initial guess
        gsl_vector *x0 = gsl_vector_alloc(2);
        gsl_vector_set(x0,0,theta0);
        gsl_vector_set(x0,1,asin(sinchi0));
        
        //Initiate solver
        const gsl_multiroot_fsolver_type *fs_t = gsl_multiroot_fsolver_dnewton;
        gsl_multiroot_fsolver *fs = gsl_multiroot_fsolver_alloc(fs_t,2);
        gsl_multiroot_function mrf = {rootfunc,2,&mi};
        gsl_multiroot_fsolver_set(fs,&mrf,x0);
        gsl_set_error_handler_off();
        
        int status, count = 0;
        
        //Iterate solver
        while(count != maxiter){
            status = gsl_multiroot_fsolver_iterate(fs);
            if(status != GSL_SUCCESS) break; //Some error has occurred
            
            gsl_vector *f = gsl_multiroot_fsolver_f(fs);
            if(gsl_multiroot_test_residual(f,rtol) != GSL_CONTINUE) break; //Done!
            
            count += 1;
        }
        
        //Free memory
        gsl_multiroot_fsolver_free(fs);
        gsl_vector_free(x0);
        
        //In case of error or aborted computation, indicate no result
        if((status != GSL_SUCCESS) || (count == maxiter)) return 1;
        
        return 0;
}

