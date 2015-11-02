# #############################################################################
# #############################################################################
#Julia wrapper for C mode solver

#This file contains the Julia wrapper for running the internal mode solver written in C (using the GSL library). The following function is provided for the mode searcher:

#   modebounces = findmode_gsl(bnd::Boundary,idx::RefractiveIndex, init::Array{Float64,1},order::Int64;rtol::Float64=1e-5,maxiter::Int64=40)
#   Use Newton's method to search for a nearby value of (<theta>,<chi>) or (angular position,angle of incidence) close to the initial trial value <init> (entered as the array [<theta0>,<chi0>]), such that the trajectory of a ray beginning in that configuration is periodic with <order> bounces (that is on the <order>th following bounce the ray returns to within a distance measure tolerance <rtol> (default 1e-5) of (<theta>,<chi>)). 
#   Returns the 2D array <modebounces> which has as its rows the (<theta>,<chi>) or (angular position,angle of incidence) values of each bounce in one cycle of a periodic mode.
#   The search is aborted and a NaN array is returned if there is a Jacobian computation error, oscillatory solution error, or <maxiter> iterations (default 40) has passed without finding a solution.
#   This function directly runs the GSL rootfinding algorithms on the ray solver C code, since the Julia wrapping only provides data recording and primitive analysis which is not required here.


# #############################################################################
# #############################################################################
#Initiate

#Dependencies
#require("basic.jl")
#require("boundary.jl")
#require("refractive_index.jl")
#require("ray_solver.jl")


# #############################################################################
# #############################################################################
#Wrapper for GSL mode finding

#Search for periodic modes (GSL multidimensional root solver)
#Use Newton's method to search for a nearby mode
function findmode_gsl(bnd::Boundary,idx::RefractiveIndex,init::Array{Float64,1},order::Int64;rtol::Float64=1e-5,maxiter::Int64=40)
    
    #Preallocate results array
    modebounces::Array{Float64,1} = Array(Float64,2*order)
    
    #Make C call
    status::Int32 = ccall((:findmode,"libcavchaos.so"), Cint,
        #Tuple of argument types
        (Clong,Ptr{Cdouble},
         Cdouble,Cdouble,
         Cdouble,Clong,
         Any,Ptr{Void},Ptr{Void},Any,Ptr{Void},Ptr{Void}),
        #Argument values
        order,modebounces, #mode-specific variables
        init[1],init[2], #initial conditions
        rtol,maxiter, #root-finding algorithm parameters
        bnd,rfunc_p,rsys!_p,idx,nfunc_p,nderiv!_p) #Cavity parameters
    
    #Return NaN on abortions or false positives (almost identical points at high chi)
    if (status == 1) || (var(modebounces[1:order]) < 1e-3)
        return nans(1,1)::Array{Float64,2}
    end
    
    #Return array where the rows are the (theta,chi) values of each bounce along the trajectory in one cycle of a periodic mode
    return reshape(modebounces,order,2)::Array{Float64,2}
    
end
