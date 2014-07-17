#require("cavbnds.jl");
require("cavidxs.jl");

# Notation:
# (r,theta) are polar position coordinates
# phi is the angle the direction of travel intersects the x-axis
# (p_r,p_theta) are "momentum" coordinates corresponding to (r,theta)


#==============================================================================
#Interface to GSL Ray Tracing


#Generate C-callable functions
function rfunc_(bndptr::Ptr{Void},theta::Cdouble)
    bnd::Boundary = unsafe_pointer_to_objref(bndptr);
    return rfunc(bnd,theta)::Cdouble;
end
function rsys!_(bndptr::Ptr{Void},theta::Cdouble,resultptr::Ptr{Cdouble})
    bnd::Boundary = unsafe_pointer_to_objref(bndptr);
    result::Array{Float64,1} = pointer_to_array(resultptr,2);
    rsys!(bnd,theta,result);
    return nothing;
end
function nfunc_(idxptr::Ptr{Void},r::Cdouble,theta::Cdouble)
    idx::RefractiveIndex = unsafe_pointer_to_objref(idxptr);
    return nfunc(idx,r,theta)::Cdouble;
end
function nsys!_(idxptr::Ptr{Void},r::Cdouble,theta::Cdouble,resultptr::Ptr{Cdouble})
    idx::RefractiveIndex = unsafe_pointer_to_objref(idxptr);
    result::Array{Float64,1} = pointer_to_array(resultptr,3);
    nsys!(idx,r,theta,result);
    return nothing;
end


#Ray traversal (GSL solver)
function rayevolve_gsl(bnd::Boundary,idx::RefractiveIndex,init::Array{Float64,1};tmax::Float64=200.0,bouncemax::Int32=int32(500),reltol::Float64=1.0e-12,abstol::Float64=1.0e-12)
    
    #Get C pointers
    const rfunc_c = cfunction(rfunc_,Cdouble,(Ptr{Void},Cdouble));
    const rsys!_c = cfunction(rsys!_,Void,(Ptr{Void},Cdouble,Ptr{Cdouble}));
    const nfunc_c = cfunction(nfunc_,Cdouble,(Ptr{Void},Cdouble,Cdouble));
    const nsys!_c = cfunction(nsys!_,Void,(Ptr{Void},Cdouble,Cdouble,Ptr{Cdouble}));
    
    #Initialize
    n::Float64 = nfunc(idx,init[1],init[2]);
    pr0::Float64 = n*cos(init[3]-init[2]);
    ptheta0::Float64 = n*init[1]*sin(init[3]-init[2]);
    
    prealloc::Int64 = int64(250*ceil(tmax));
    raypath1::Array{Float64,1} = Array(Float64,prealloc);
    raypath2::Array{Float64,1} = Array(Float64,prealloc);
    bounces1::Array{Int32,1} = Array(Int32,bouncemax);
    bounces2::Array{Float64,1} = Array(Float64,bouncemax);
    lengths::Array{Int32,1} = Array(Int32,2);
    
    #Make C call
    ccall((:rayevolve,"./libray.so"), Void,
            
            (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cint},Ptr{Cdouble},Ptr{Cint},
             Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cint,Cdouble,Cdouble,
             Any,Ptr{Void},Ptr{Void},Any,Ptr{Void},Ptr{Void}),
            
            raypath1,raypath2,bounces1,bounces2,lengths,
            init[1],init[2],pr0,ptheta0,tmax,bouncemax,reltol,abstol,
            bnd,rfunc_c,rsys!_c,idx,nfunc_c,nsys!_c);
    
    #Truncate allocated arrays for output
    resize!(raypath1,lengths[1]); resize!(raypath2,lengths[1]);
    resize!(bounces1,lengths[2]); resize!(bounces2,lengths[2]);
    
    return ([raypath1 raypath2]::Array{Float64,2},bounces1::Array{Int32,1},[raypath2[bounces1] bounces2]::Array{Float64,2});
    
end
