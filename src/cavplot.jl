#require("cavbnds.jl");
require("cavidxs.jl");
import PyPlot.plt;
using PyCall;

#Add current directory to python import path
unshift!(PyVector(pyimport("sys")["path"]), "");
interact = pyimport("interact");

mpl = PyPlot.matplotlib;
plt.ioff();

#==============================================================================
#Plotting tools


#Plot the boundary
function plotbnd(bnd::Boundary;axes=plt.gca(),normvec::Bool=false)
    #Plots a boundary, with the option of showing its normal vectors as red emanating lines.
    theta::Array{Float64,1} = linspace(0,2*pi,500);
    r::Array{Float64,1} = bnd.rfunc(theta);
    plt.plot(r.*cos(theta),r.*sin(theta),"g",axes=axes);
    plt.axis([-r[1],r[1],-r[1],r[1]],axes=axes);
    
    if normvec
        #Plot the norm vector along boundary
        theta_n::Array{Float64,1} = linspace(0,2*pi,51);
        r_n::Array{Float64,1},normang::Array{Float64,1} = bnd.sys(theta_n);
        for i = 1:length(theta_n)
            plt.plot([r_n[i]*cos(theta_n[i]),r_n[i]*cos(theta_n[i])+0.1*cos(normang[i])],[r_n[i]*sin(theta_n[i]),r_n[i]*sin(theta_n[i])+0.1*sin(normang[i])],"r",axes=axes);
        end
    end
    
    return nothing;
end


#Plot the refractive index distribution
function plotidx(idx::RefractiveIndex,bnd::Boundary;axes=plt.gca())
    #Plots a filled labelled contour plot of the cavity interior's refractive index distribution, where darker shades indicate higher refractive index.
    thetares::Int64 = 100; rres::Int64 = 100;
    
    theta::Array{Float64,2} = transpose(linspace(0,2*pi,thetares));
    rbnd::Array{Float64,2} = bnd.rfunc(theta);
    rs::Array{Float64,2} = linspace(0.0,1.0,rres)*rbnd;
    thetas::Array{Float64,2} = ones(Float64,rres)*theta;
    
    x::Array{Float64,2} = rs.*cos(thetas);
    y::Array{Float64,2} = rs.*sin(thetas);
    n::Array{Float64,2} = idx.nfunc(rs,thetas);
    
    plt.contourf(x,y,n,100,cmap="bone_r",axes=axes);
    contobj = plt.contour(x,y,n,100,cmap="bone_r",alpha=0.2,axes=axes);
    levels = contobj[:levels][1:10:end];
    plt.clabel(contobj,levels,colors="k");
    return nothing;
end
function plotidx(idx::RefractiveIndex;axes=plt.gca())
    #Plots a filled labelled contour plot of the refractive index distribution
    res::Int64 = 200;
    x::Array{Float64,2} = ones(Float64,res)*transpose(linspace(-1.0,1.0,res));
    y::Array{Float64,2} = linspace(-1.0,1.0,res)*ones(Float64,(1,res));
    rs::Array{Float64,2} = hypot(x,y);
    thetas::Array{Float64,2} = atan2(y,x);
    n::Array{Float64,2} = idx.nfunc(rs,thetas);
    
    plt.contourf(x,y,n,100,cmap="bone_r",axes=axes);
    contobj = plt.contour(x,y,n,100,cmap="bone_r",alpha=0.2,axes=axes);
    levels = contobj[:levels][1:10:end];
    plt.clabel(contobj,levels,colors="k");
    plt.title("Plot of $idx");
    return nothing;
end


#Draw rays
function plotrays(raypath::Array{Float64,2};axes=plt.gca())
    plt.plot(raypath[:,1].*cos(raypath[:,2]),raypath[:,1].*sin(raypath[:,2]),"y-",lw=0.4,aa=false,axes=axes);
    return nothing;
end


#Plot Poincare Surface of Section
function plotpss(bouncesval::Array{Array{Float64,2},1},bnd::Boundary=nothing,idx::RefractiveIndex=nothing;axes=plt.gca())
    mpl[:rc]("text", usetex=true);

    #Emphasize first set of points
    bounces::Array{Float64,2} = bouncesval[1];
    plt.plot(bounces[:,1],sin(bounces[:,2]),".",ms=8,mec="k",mfc="r",mew=1.5,picker=true,zorder=length(bouncesval),axes=axes);

    #Scatter remaining points
    for i=2:length(bouncesval)
        bounces = bouncesval[i];
        uclst::Float64,order::Int64 = interact[:unclusteredness](bounces);
        color = interact[:colorfunc](uclst);
        plt.plot(bounces[:,1],sin(bounces[:,2]),".",color=color,markersize=2,picker=true,axes=axes);
    end
    plt.axis([-pi,pi,-1,1]);
    plt.xlabel("\$\\displaystyle\\theta\$", fontsize=18);
    plt.ylabel("\$\\displaystyle\\sin\\chi\$", fontsize=18);

    #Show critical angle
    if bnd != nothing && idx != nothing
        theta = linspace(-pi,pi,200);
        sinchi0 = 1.0/idx.nfunc(float64(bnd.rfunc(theta)),theta);
        plt.plot(theta,sinchi0,"m",theta,-sinchi0,"m",axes=axes);
    end
end

