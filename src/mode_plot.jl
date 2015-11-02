# #############################################################################
# #############################################################################
#Mode Plotting functions

#This file contains functions for plotting modes in a cavity. 

# #############################################################################
# #############################################################################
#Initiate

#require("util.jl")
#require("boundary.jl")
#require("refractive_index.jl")
#require("io.jl")
#require("ray_solver.jl")
#require("ray_plot.jl")
#require("clustering.jl")


# #############################################################################
# #############################################################################
#Preallocated array based non-interactive PSS plot

#Get corresponding index of a position in an array representing regions in a range
#This is a nearest neighbor interpolation for points within the range
function getbinindex(pos::Float64,min::Float64,max::Float64,res::Int64)
    #Check if outside region
    if (pos < min && pos < max) || (pos > min && pos > max)
        return 0
    else
        #Return the bin number containing pos, if there are res bins between
        #min and max
        return convert(Int64,cld(res*(pos-min),max-min))
    end
end

#Incorporate bouncepoint dataset into arrays, letting each point replace those in
#its corresponding bin if it has a higher priority
function addbouncepointstopss(bouncepoints::Array{Float64,2},resultid::Int64,pssarray::Array{Float64,2},residarray::Array{Int64,2},range::Array{Float64,1})
    #Region parameters
    #Range is [thetamin,thetamax,chimin,chimax]
    #assert(size(pssarray) == size(residarray))

    #Determine score for set of bouncepoints based on degree of cluster spread
    #Each nonnegative cluster spread value is scored with a value in (0.0,1.0), stationary at small
    #distances, and asymptotic to 1 as it approaches infinity. Lower score has higher priority.
    #Cluster spread estimated using variance measures consistent with metric on pss
    order::Int64,clustsize::Float64 = findcluster(bouncepoints,maxtestinterval=25)
    clusterstddev::Float64 = sqrt(var(sin(bouncepoints[1:order:end,1]))
                                  + var(cos(bouncepoints[1:order:end,1]))
                                  + var(map(chi2whisporder,bouncepoints[1:order:end,2])))
    score::Float64 = sq(1.0 - 1.0 / (1.0 + 10*clusterstddev))
    
    #Loop over all bouncepoints
    for i = 1:size(bouncepoints,1)
        #Get corresponding index in array
        thetaindex::Int64 = getbinindex(bouncepoints[i,1],range[1],range[2],size(pssarray,2))
        chiindex::Int64 = getbinindex(bouncepoints[i,2],range[4],range[3],size(pssarray,1))
        
        #Skip if out of range
        if thetaindex == 0 || chiindex == 0
            continue
        end
        
        #Replace value if score is lower
        if score < pssarray[chiindex,thetaindex]
            pssarray[chiindex,thetaindex] = score
            residarray[chiindex,thetaindex] = resultid
        end
    end
end

#Incorporate all bouncepoint datasets in a directory
function adddirtopss(resultsdir::AbstractString;range::Array{Float64,1}=[-pi,pi,-pi/6,pi/2],res::Int64=400)
    #Initialize Arrays
    pssarray::Array{Float64,2} = fill(10.0,res,res) #array of some large number by default
    residarray::Array{Int64,2} = zeros(Int64,res,res)

    #Incorporate data
    resultids::Array{Int64,1} = collectresultids(resultsdir,"bouncepoints")
    for resultid in resultids
        bouncepoints::Array{Float64,2} = readresults(resultsdir,resultid,"bouncepoints")
        addbouncepointstopss(bouncepoints,resultid,pssarray,residarray,range)
    end

    #Return array of pss scores, and array of result ids
    return pssarray,residarray
end

#Plot PSS
function plotpss(resultsdir::AbstractString;axes=plt[:gca](),range::Array{Float64,1}=[-pi,pi,-pi/6,pi/2],res::Int64=400,colorrange::Array{Float64,1}=[-0.1,1.1])
    #Colormap normalization
    colornorm::PyObject = plt[:matplotlib][:colors][:Normalize](colorrange[1],colorrange[2])
    
    #Get PSS
    pssarray::Array{Float64,2},residarray::Array{Int64,2} = adddirtopss(resultsdir;range=range,res=res)
    
    #Plot pssarray
    plt[:imshow](pssarray,norm=colornorm,cmap="spectral_r",extent=range,aspect=3)
    plt[:title]("Poincare Surface of Section")
    plt[:xlabel]("Polar Angle \$\\theta\$")
    plt[:ylabel]("Angle of Incidence \$\\chi\$")

    return pssarray,residarray
end


# #############################################################################
# #############################################################################
#PSS plot I/O

#Get PSS plot image file location
function getpssimgfile(resultsdir::AbstractString)
    return joinpath(resultsdir,"pssplot.png")
end

#Write image of PSS plot
function writepssimg(resultsdir::AbstractString;kwargs...)
    #Create new figure
    pssplot = plt[:figure](figsize=(9,9))[:gca](aspect="equal")

    #Plot PSS non-interactively
    plotpss(resultsdir;axes=pssplot,kwargs...)

    #Save image
    plt[:savefig](getpssimgfile(resultsdir))
    return nothing
end
