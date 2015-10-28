# #############################################################################
# #############################################################################
# Plotting commands

#This file contains the following functions for common plots, used to render low-res images of ray paths for easy qualitative identification of modes.

#   plotbnd(bnd::Boundary;axes=plt.gca(),normvec::Bool=false)
#   Plots the boundary <bnd> on the specified matplotlib axes <axes>. If <normvec> is true, also plots some arrows showing the normal vector at some points along the boundary, for a visual verification of the derivative functions.

#   plotidx(idx::RefractiveIndex,[bnd::Boundary];axes=plt.gca())
#   Plots a filled labelled contour plot of the refractive index distribution <idx> on the specified matplotlib axes <axes>. If the boundary <bnd> is specified, the contour plot is not plotted in the entire axes but limited to the region within the boundary, so this can be plotted alongside the boundary plot for a nice image of the cavity.

#   rayline = plotrays(raypath::Array{Float64,2};axes=plt.gca())
#   Plots on the specified matplotlib axes <axes> the rays along the path specified by <raypath>, which has rows corresponding to the polar (<r>,<theta>) coordinates of the ray position. Returns the matplotlib line object for the plotted rays.

#The following functions are used for storing or retrieving the cavity and ray plots as bitmap images, which are sufficient for obtaining a qualitative idea of any mode present, but much less disk space intensive than storing the entire raypath data directly.

#   cavityimgpath = getcavityimgfile(resultsdir::AbstractString)
#   Returns the path <cavityimgpath> to the cavity image file for a specified simulation run in the results directory <resultsdir>.

#   rayimgpath = getrayimgfile(resultsdir::AbstractString,resultid::Int64)
#   Returns the path <rayimgpath> to the ray image file for the <resultid>th computation of a specified simulation run in the results directory <resultsdir>.

#   range = writecavityimg(resultsdir::AbstractString,bnd::Boundary,idx::RefractiveIndex)
#   Writes an image of the cavity with boundary <bnd> and index distribution <idx> in the results directory <resultsdir>. Returns the range <range> of the matplotlib axes plot for calibrating the ray images.

#   rayline = writerayimg(resultsdir::AbstractString,resultid::Int64, raypath::Array{Float64,2},range::Array{Float64,1}=[-1.0,1.0,-1.0,1.0])
#   Writes an image of the rays following the path <raypath> from the <resultid>th computation of a specified simulation run in the results directory <resultsdir>, adjusting the matplotlib axes range to <range>. Returns the matplotlib line object <rayline> from the plot.

#   writerayimg(resultsdir::AbstractString,resultid::Int64,raypath::Array{Float64,2},rayline::PyObject)
#   Writes an image of the rays following the path <raypath> from the <resultid>th computation of a specified simulation run in the results directory <resultsdir>. Obtains image by altering the line data in the matplotlib line object <rayline> of an existing ray image, for a speedier plot.


# #############################################################################
# #############################################################################
#Initiate

#Dependencies
#require("util.jl")
#require("boundary.jl")
#require("refractive_index.jl")
#require("io.jl")
import PyPlot.plt
plt.ioff()
using PyCall


# #############################################################################
# #############################################################################
#Plotting tools

#Common preparation of empty plot for cavity and ray plotting functions
function prepareplot(;figwidth::Float64=5.0,figheight::Float64=5.0)
    ax = plt.figure(figsize=(figwidth,figheight))[:gca](aspect="equal")
    plt.gcf()[:tight_layout]()
    return ax
end

#Plot the boundary
function plotbnd(bnd::Boundary;axes=plt.gca(),normvec::Bool=false)
    #Plots a boundary, with the option of showing its normal vectors as red emanating lines.
    theta::Array{Float64,1} = linspace(0,2*pi,500)
    r::Array{Float64,1} = rfunc(bnd,theta)
    plt.plot(r.*cos(theta),r.*sin(theta),"g",axes=axes)
    plt.axis([-r[1],r[1],-r[1],r[1]],axes=axes)
    
    if normvec
        #Plot the norm vector along boundary
        theta_n::Array{Float64,1} = linspace(0,2*pi,51)
        r_n::Array{Float64,1},alpha::Array{Float64,1} = rsys(bnd,theta_n)
        for i = 1:length(theta_n)
            plt.arrow(r_n[i]*cos(theta_n[i]),r_n[i]*sin(theta_n[i]),0.1*cos(alpha[i]),0.1*sin(alpha[i]);color="r",axes=axes)
        end
    end
    
    return nothing
end


#Plot the refractive index distribution
function plotidx(idx::RefractiveIndex,bnd::Boundary;axes=plt.gca())
    #Plots a filled labelled contour plot of the cavity interior's refractive index distribution, where darker shades indicate higher refractive index.
    thetares::Int64 = 100; rres::Int64 = 100
    
    theta::Array{Float64,2} = transpose(linspace(0,2*pi,thetares))
    R::Array{Float64,2} = rfunc(bnd,theta)
    rs::Array{Float64,2} = linspace(0.0,1.0,rres)*R
    thetas::Array{Float64,2} = ones(Float64,rres)*theta
    
    x::Array{Float64,2} = rs.*cos(thetas)
    y::Array{Float64,2} = rs.*sin(thetas)
    n::Array{Float64,2} = nfunc(idx,rs,thetas)
    
    plt.contourf(x,y,n,100,cmap="bone_r",axes=axes)
    contobj = plt.contour(x,y,n,100,cmap="bone_r",alpha=0.2,axes=axes)
    levels = contobj[:levels][1:10:end]
    plt.clabel(contobj,levels,colors="k")
    return nothing
end

function plotidx(idx::RefractiveIndex;axes=plt.gca())
    #Plots a filled labelled contour plot of the refractive index distribution
    res::Int64 = 200
    x::Array{Float64,2} = ones(Float64,res)*transpose(linspace(-1.0,1.0,res))
    y::Array{Float64,2} = linspace(-1.0,1.0,res)*ones(Float64,(1,res))
    rs::Array{Float64,2} = hypot(x,y)
    thetas::Array{Float64,2} = atan2(y,x)
    n::Array{Float64,2} = nfunc(idx,rs,thetas)
    
    plt.contourf(x,y,n,100,cmap="bone_r",axes=axes)
    contobj = plt.contour(x,y,n,100,cmap="bone_r",alpha=0.2,axes=axes)
    levels = contobj[:levels][1:10:end]
    plt.clabel(contobj,levels,colors="k")
    plt.title("Plot of $idx")
    return nothing
end


#Draw rays
function plotrays(raypath::Array{Float64,2};axes=plt.gca())
    return plt.plot(raypath[:,1].*cos(raypath[:,2]),raypath[:,1].*sin(raypath[:,2]),"y-",lw=0.4,aa=false,axes=axes)[1]
end


# #############################################################################
# #############################################################################
#Ray image I/O

#Get cavity image file location
function getcavityimgfile(resultsdir::AbstractString)
    return joinpath(resultsdir,"cavity.png")
end

#Get ray image file location
function getrayimgfile(resultsdir::AbstractString,resultid::Int64)
    return joinpath(resultsdir,getresultfname(resultid,"ray")*".png")
end

#Attempts to store image of cavity if not already present. Returns the axis range of the cavity.
function writecavityimg(resultsdir::AbstractString,bnd::Boundary,idx::RefractiveIndex)
    #Create new figure
    cavityplot = prepareplot()
    
    #Plot boundary to obtain axis range
    plotbnd(bnd,axes=cavityplot)
    range = [plt.axis()...]
    
    #Continue plotting if cavity image file does not already exist
    cavityimgfile = getcavityimgfile(resultsdir)
    if !isfile(cavityimgfile)
        #Plot index distribution
        plotidx(idx,bnd,axes=cavityplot)
        plt.savefig(cavityimgfile)
    end
    
    plt.close()
    return range
end

#Stores newly created image of ray and returns the matplotlib line object of the ray plot for more efficient future editing.
function writerayimg(resultsdir::AbstractString,resultid::Int64,raypath::Array{Float64,2},range::Array{Float64,1}=[-1.0,1.0,-1.0,1.0])
    #Create new figure
    rayplot = prepareplot()
    
    #Plot
    rayline = plotrays(raypath,axes=rayplot)
    plt.axis(range)
    plt.axis("off")
    
    #Save image
    plt.savefig(getrayimgfile(resultsdir,resultid),transparent=true)
    return rayline
end

#Stores image of ray after editing the matplotlib line object of an existing ray plot.
function writerayimg(resultsdir::AbstractString,resultid::Int64,raypath::Array{Float64,2},rayline::PyObject)
    #Edit existing ray line
    rayline[:set_xdata](raypath[:,1].*cos(raypath[:,2]))
    rayline[:set_ydata](raypath[:,1].*sin(raypath[:,2]))
    
    #Save image
    plt.savefig(getrayimgfile(resultsdir,resultid),transparent=true)
    return nothing
end

