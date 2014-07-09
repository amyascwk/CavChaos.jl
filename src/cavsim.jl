#require("cavbnds.jl");
require("cavidxs.jl");

# Notation:
# (r,theta) are polar position coordinates
# phi is the angle the direction of travel intersects the x-axis
# (p_r,p_theta) are "momentum" coordinates corresponding to (r,theta)


#==============================================================================
#Native Julia Ray Tracing


#Compute "momentum" p coordinates from direction angle phi
function phi2p(r::Float64,theta::Float64,phi::Float64,idx::RefractiveIndex)
    #p_r,p_theta obtained from eikonal eqn (See Hakan p.41)
    #(\del S)^2 = n^2 = p_r^2 + \frac{p_\theta^2}{r^2}
    n::Float64 = idx.nfunc(r,theta);
    p_r::Float64 = n*cos(phi-theta);
    p_theta::Float64 = n*r*sin(phi-theta);
    return (p_r,p_theta);
end

# ******** Warning: Julia solver is obsolete ********

# #Ray reflection (Julia solver)
# function raybounce(r0::Float64,theta0::Float64,r::Float64,theta::Float64,bnd::Boundary,idx::RefractiveIndex)
#     #Re-compute a bounced ray coordinate vector
#     # r0, theta0:               Pre-collision position
#     # r,theta,p_r,p_theta:      Post-collision coordinate vector
# 
    #Binary search for intersection of trajectory with cavity boundary
#     x0::Float64 = r0*cos(theta0); y0::Float64 = r0*sin(theta0);
#     xA::Float64 = x0; yA::Float64 = y0;
#     xB::Float64 = r*cos(theta); yB::Float64 = r*sin(theta);
#     while true
#         xC::Float64 = 0.5*(xA+xB); yC::Float64 = 0.5*(yA+yB);
#         rC::Float64 = hypot(xC,yC); thetaC::Float64 = atan2(yC,xC);
#         RC::Float64 = bnd.rfunc(thetaC);
#         if abs(rC-RC) < 1e-12
#             break;
#         elseif rC > RC
#             xB = xC; yB = yC;
#         else
#             xA = xC; yA = yC;
#         end
#     end
# 
    #Compute reflected ray
#     RC::Float64,normang::Float64 = bnd.sys(thetac);
#     chi::Float64 = atan2(yC-y0,xC-x0) - normang; #angle of incidence
#     phi::Float64 = pi - chi + normang; #angle of reflected light to horizontal
#     prc::Float64,pthetac::Float64 = phi2p(rC,thetaC,phi,idx);
#     
#     S::Array{Float64,1} = [rC,thetaC,prc,pthetac];
#     bounce::Array{Float64,2} = [thetaC chi];
#     
#     return (S,bounce);
# end
# 

# #Ray traversal through medium
# function rayevolve(r0::Float64=1.0,theta0::Float64=0.0,phi0::Float64=0.5;bnd::Boundary=CircularBnd(),idx::RefractiveIndex=UniformIdx(),dt::Float64=0.001,tmax::Float64=10.0,reltol::Float64 = 1.0e-12,abstol::Float64 = 1.0e-12)
#     #Numerically integrate the ray by solving the Eikonal equation
#     #r0, theta0, ptheta0:    Initial values of r, theta, and ptheta.
#     #tmax:                   Total time for integration.
#     #dt:                     Time step for the integration.
#     
#     
#     function F(t::Float64,y::Array{Float64,1})
#         #The derivative function f for the traversal of the ray, such that
#         #    y' = f(y)
#         #where y = [r,theta,pr,ptheta] is 
#         #expresses the parameters of the ray at a given point.
#         #(See (2.15) and (2.16) in Hakan's thesis)
#         #    r' = \frac{p_r}{n}, \theta' = \frac{p_\theta}{r^2 n},
#         #    p_r' = \frac{\partial n}{\partial r} + \frac{p_\theta^2}{n r^3},
#         #    p_\theta' = \frac{\partial n}{\partial \theta}
#         n::Float64,dr_n::Float64,dtheta_n::Float64 = idx.sys(y[1],y[2]);
#         dydt::Array{Float64,1} = [y[3]/n, y[4]/(n*sq(y[1])), dr_n+sq(y[4])/(n*cb(y[1])), dtheta_n];
#         return dydt;
#     end
#     
#     #Initiation
#     p_r0::Float64,p_theta0::Float64 = phi2p(r0,theta0,phi0,idx);
#     S0::Array{Float64,1} = [r0, theta0, p_r0, p_theta0];
#     
#     #Adapted from ODE.jl ------------------------------------------------------
#     # ************* IMPROVEMENT: Replace with GSL for speed? ******************
#     #Dormand-Prince solver
#     const a::Array{Float64,2} = [ 0 0 0 0 0 0;
#           1/5 0 0 0 0 0;
#           3/40 9/40 0 0 0 0;
#           44/45 -56/15 32/9 0 0 0;
#           19372/6561 -25360/2187 64448/6561 -212/729 0 0;
#           9017/3168 -355/33 46732/5247 49/176 -5103/18656 0;
#           35/384 0 500/1113 125/192 -2187/6784 11/84];
#     const bs::Array{Float64,1} = [5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40];
#     const bp::Array{Float64,1} = [35/384,0,500/1113,125/192,-2187/6784,11/84,0];
#     
#     const c::Array{Float64,2} = sum(a, 2); # consistency condition
# 
    # Initialization
#     h::Float64 = dt; # initial guess at a step size
#     hmax::Float64 = 10*dt;
#     hmin::Float64 = 1e-5*dt;
#     t::Float64 = 0.0;
#     S::Array{Float64,1} = S0;
#     rval::Array{Float64,1} = [S[1]];
#     thetaval::Array{Float64,1} = [S[2]];
#     bounces::Array{Float64,2} = Array(Float64,200,2);
#     bouncenum::Int64 = 0;
# 
#     k::Array{Float64,2} = Array(Float64,(4,length(c)));
#     k[:,1] = F(t,S);
# 
#     while t < tmax && h >= hmin
# 
        #(p-1)th and pth order estimates
#         Ss::Array{Float64,1} = S + h*bs[1]*k[:,1];
#         Sp::Array{Float64,1} = S + h*bp[1]*k[:,1];
#         for j = 2:length(c)
#             dS::Array{Float64,1} = a[j,1]*k[:,1];
#             for i = 2:j-1
#                 dS += a[j,i]*k[:,i];
#             end
#             k[:,j] = F(t + h*c[j], S + h*dS);
# 
            # compute the (p-1)th order estimate
#             Ss = Ss + h*bs[j]*k[:,j];
#             # compute the pth order estimate
#             Sp = Sp + h*bp[j]*k[:,j];
#         end
# 
        # estimate the local truncation error
#         gamma1::Array{Float64,1} = Ss - Sp;
# 
        # Estimate the error and the acceptable error
#         delta::Float64 = max(norm(gamma1, Inf)) # actual error
#         tau::Float64 = max(reltol*norm(S,Inf),abstol) # allowable error
#         
#         #Check for acceptability of error
#         if delta <= tau
#             
#             #Verify conservation of Hamiltonian function
#             if mod(length(rval),1000) == 0
#                 run(`clear`);
#                 println("t = ",t);
#                 H::Float64 = sq(S[3])+sq(S[4])/sq(S[1]) - sq(idx.nfunc(S[1],S[2]));
#                 if abs(H) > 1e-9
#                     println("Warning: Hamiltonian is ",H,", might not be conserved.");
#                 end
#             end
#         
#             #Check for boundary crossing
#             dr::Float64 = Sp[1] - bnd.rfunc(Sp[2]); #Boundary crossing error
#             if dr > 0
#                 Sp, bounces[bouncenum+1,:] = raybounce(S[1:2]...,Sp[1:2]...,bnd,idx);
#                 bouncenum += 1;
#             end
#             
#             # Update the solution
#             t = t + h;
#             S = Sp;
#             push!(rval,S[1]);
#             push!(thetaval,S[2]);
#             
#             if c[end] == 1
#                 k[:,1] = k[:,end];
#             else
#                 k[:,1] = F(t,S); # first stage
#             end
#         end
# 
        # Update the step size
#         h = min(hmax, 0.8*h*(tau/delta)^0.2);
#         
#     end
#     
#     if t < tmax
#       println("Step size grew too small. t=", t, ", h=", h, ", S=", S)
#     end
#     # -------------------------------------------------------------------------
# 
#     return (bounces[1:bouncenum,:],rval,thetaval)
# end


#==============================================================================
#Interface to GSL Ray Tracing

#Make functions top level in order to be C-callable
function rayev_1arg(func_p::Ptr{Void},theta::Cdouble)
    func = unsafe_pointer_to_objref(func_p)::Function;
    return func(theta)::Cdouble;
end
function rayev_1arg_(func_p::Ptr{Void},theta::Cdouble,array_p::Ptr{Cdouble})
    func! = unsafe_pointer_to_objref(func_p)::Function;
    arr = pointer_to_array(array_p,3);
    func!(theta,arr);
    return nothing;
end
function rayev_2arg(func_p::Ptr{Void},r::Cdouble,theta::Cdouble)
    func = unsafe_pointer_to_objref(func_p)::Function;
    return func(r,theta)::Cdouble;
end
function rayev_2arg_(func_p::Ptr{Void},r::Cdouble,theta::Cdouble,array_p::Ptr{Cdouble})
    func! = unsafe_pointer_to_objref(func_p)::Function;
    arr = pointer_to_array(array_p,3);
    func!(r,theta,arr);
    return nothing;
end


#Ray traversal (GSL solver)
function rayevolve_gsl(r0::Float64=1.0,theta0::Float64=0.0,phi0::Float64=0.5;bnd::Boundary=CircularBnd(),idx::RefractiveIndex=UniformIdx(),tmax::Float64=200.0,reltol::Float64 = 1.0e-12,abstol::Float64 = 1.0e-12)
    
    #Get C pointers
    const rfunc_c = cfunction(rayev_1arg,Cdouble,(Ptr{Void},Cdouble));
    const rsys_c = cfunction(rayev_1arg_,Void,(Ptr{Void},Cdouble,Ptr{Cdouble}));
    const nfunc_c = cfunction(rayev_2arg,Cdouble,(Ptr{Void},Cdouble,Cdouble));
    const nsys_c = cfunction(rayev_2arg_,Void,(Ptr{Void},Cdouble,Cdouble,Ptr{Cdouble}));
    
    #Initialize
    pr0::Float64, ptheta0::Float64 = phi2p(r0,theta0,phi0,idx);
    prealloc::Int64 = int64(250*ceil(tmax));
    raypath1::Array{Float64,1} = Array(Float64,prealloc);
    raypath2::Array{Float64,1} = Array(Float64,prealloc);
    bounces1::Array{Float64,1} = Array(Float64,prealloc);
    bounces2::Array{Float64,1} = Array(Float64,prealloc);
    
    #Make C call
    ccall((:rayevolve,"./libray.so"), Void,
    
            (Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},
             Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,
             Any,Ptr{Void},Any,Ptr{Void},Any,Ptr{Void},Any,Ptr{Void}),
            
            raypath1,raypath2,bounces1,bounces2,
            r0,theta0,pr0,ptheta0,tmax,reltol,abstol,
            bnd.rfunc,rfunc_c,bnd.sys!,rsys_c,idx.nfunc,nfunc_c,idx.sys!,nsys_c);
    
    #Truncate output
    stepnum::Int64 = int64(raypath1[1]-1);
    bouncenum::Int64 = int64(bounces1[1]-1);
    return ([raypath1[2:stepnum] raypath2[2:stepnum]],[bounces1[2:bouncenum] bounces2[2:bouncenum]]);
    
end
