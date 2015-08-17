# #############################################################################
# #############################################################################
#Julia wrappers for C solvers, and ray solver

#This file contains the Julia wrappers and functions for running the internal solvers written in C (using the GSL library). The following Julia function is provided for the ray solver:

#   raypath,bounce_indices,bounce_points = rayevolve_gsl(bnd::Boundary, idx::RefractiveIndex,init::Array{Float64,1};tmax::Float64=200.0,bouncemax::Int64=500,reltol::Float64=1.0e-12,abstol::Float64=1.0e-12)
#   Computes the trajectory and bounce points in a cavity with boundary <bnd> and index distribution <idx>, of a ray with initial conditions <init>, which is an array of initial (<r>,<theta>,<phi>) or (radial position, angular position, angle of propagation) values or the ray. Runs the ODE solver with relative and absolute tolerances <reltol> and <abstol> (default 1.0e-12 each), and halts computation when pathlength reaches <tmax> (default 200.0) or number of bounces reaches <bouncemax> (default 500).
#   Returns the 2D array <raypath> of (<r>,<theta>) positions traversed in the ray's path as its rows, the 1D array <bounces_indices> of row indices in <raypath> that bounces occur, and the 2D array <bounce_points> with the (<theta>,<chi>) or (angular position,angle of incidence) values of each bounce as its rows.
#   This means that <raypath>[<bounces_indices>,2] == <bounce_points>[:,1].


# #############################################################################
# #############################################################################
#Initiate

#Dependencies
#require("basic.jl")
#require("boundary.jl")
#require("refractive_index.jl")


# #############################################################################
# #############################################################################
#Interface for Julia callbacks in C

#Generate C-callable functions wrapping the relevant julia functions for boundary and refractive index types
function rfunc_c(bndptr::Ptr{Void},theta::Cdouble)
    bnd::Boundary = unsafe_pointer_to_objref(bndptr)
    return rfunc(bnd,theta)::Cdouble
end
function rsys!_c(bndptr::Ptr{Void},theta::Cdouble,resultptr::Ptr{Cdouble})
    bnd::Boundary = unsafe_pointer_to_objref(bndptr)
    result::Array{Float64,1} = pointer_to_array(resultptr,2)
    rsys!(bnd,theta,result)
    return nothing
end
function nfunc_c(idxptr::Ptr{Void},r::Cdouble,theta::Cdouble)
    idx::RefractiveIndex = unsafe_pointer_to_objref(idxptr)
    return nfunc(idx,r,theta)::Cdouble
end
function nderiv!_c(idxptr::Ptr{Void},r::Cdouble,theta::Cdouble,resultptr::Ptr{Cdouble})
    idx::RefractiveIndex = unsafe_pointer_to_objref(idxptr)
    result::Array{Float64,1} = pointer_to_array(resultptr,3)
    nderiv!(idx,r,theta,result)
    return nothing
end

#Get C function pointers
const rfunc_p = cfunction(rfunc_c,Cdouble,(Ptr{Void},Cdouble))
const rsys!_p = cfunction(rsys!_c,Void,(Ptr{Void},Cdouble,Ptr{Cdouble}))
const nfunc_p = cfunction(nfunc_c,Cdouble,(Ptr{Void},Cdouble,Cdouble))
const nderiv!_p = cfunction(nderiv!_c,Void,(Ptr{Void},Cdouble,Cdouble,Ptr{Cdouble}))


# #############################################################################
# #############################################################################
#Wrapper for GSL Ray Tracing

#Ray traversal (GSL ODE solver)
#Computes the trajectory and bounce points of a ray
function rayevolve_gsl(bnd::Boundary,idx::RefractiveIndex,init::Array{Float64,1};tmax::Float64=200.0,bouncemax::Int64=500,reltol::Float64=1.0e-12,abstol::Float64=1.0e-12)
    
    #Convert initial (r,theta,phi) conditions to the vector (r,theta,pr,ptheta)
    # Notation:
    # (r,theta) are polar position coordinates
    # phi is the angle the direction of travel makes with the horizontal
    # (pr,ptheta) are "momentum" coordinates corresponding to (r,theta)
    # chi is the angle of incidence at each bounce
    const n::Float64 = nfunc(idx,init[1],init[2])
    const pr0::Float64 = n*cos(init[3]-init[2])
    const ptheta0::Float64 = n*init[1]*sin(init[3]-init[2])
    
    #Preallocate results arrays -----------------------------------------------
    #Estimate length required
    const prealloc::Int64 = int64(250*ceil(tmax)) #Length to preallocate raypath array
    #Initialize arrays for the 1st and 2nd columns of raypath matrix storing 
    #(r,theta) positions along the path of the ray as rows
    raypath_r::Array{Float64,1} = Array(Float64,prealloc) #raypath r value
    raypath_theta::Array{Float64,1} = Array(Float64,prealloc) #raypath theta value
    #Initialize arrays for the 1st and 2nd columns of the bounce information matrix 
    #storing the corresponding row of the raypath array and angle of incidence (chi) 
    #values in each row
    bounce_indices::Array{Int32,1} = Array(Int32,bouncemax) #bounce indices
    bounces_chi::Array{Float64,1} = Array(Float64,bouncemax) #bounce chi value
    #Record actual number of rows of the raypath and bounces arrays here
    lengths::Array{Int32,1} = Array(Int32,2)
    
    #Make C call --------------------------------------------------------------
    ccall((:rayevolve,"../lib/libcavchaos.so"), Void,
        #Tuple of argument types
        (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cint},Ptr{Cdouble},Ptr{Cint},
         Cdouble,Cdouble,Cdouble,Cdouble,
         Cdouble,Cint,Cdouble,Cdouble,
         Any,Ptr{Void},Ptr{Void},Any,Ptr{Void},Ptr{Void}),
        #Argument values
        raypath_r,raypath_theta,bounce_indices,bounces_chi,lengths, #Storage arrays
        init[1],init[2],pr0,ptheta0, #Initial conditions
        tmax,int32(bouncemax),reltol,abstol, #Simulation parameters
        bnd,rfunc_p,rsys!_p,idx,nfunc_p,nderiv!_p) #Cavity parameters
    
    #Truncate allocated arrays for output
    resize!(raypath_r,lengths[1])
    resize!(raypath_theta,lengths[1])
    resize!(bounce_indices,lengths[2])
    resize!(bounces_chi,lengths[2])
    
    #Return (raypath, bounce_indices, bounce_points), where the rows of raypath are the (r,theta) positions along the ray's path, bounce_indices contains the indices of raypath that bounces occur, and the rows of bounce_points are the (theta,chi) values of each bounce
    return ([raypath_r raypath_theta]::Array{Float64,2}, #raypath position array
        bounce_indices::Array{Int32,1}, #indices in raypath corresponding to bounces
        [raypath_theta[bounce_indices] bounces_chi]::Array{Float64,2}) #bounce points
    
end
