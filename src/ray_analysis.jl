# #############################################################################
# #############################################################################
#Analysis of ray

#This file contains functions for analyzing and computing properties of a given ray.

# #############################################################################
# #############################################################################
#Initiate

#require("util.jl")
#require("boundary.jl")
#require("refractive_index.jl")


# #############################################################################
# #############################################################################
#Computation of mode information

#Path integration
#Calculates the integrals \int{n(s)ds} and \int{ds} over a path, which is used to obtain the action and the path length.
function getpathintegrals(raypath::Array{Float64,2},idx::RefractiveIndex)
    #Initiate
    pathlength::Float64 = 0
    action::Float64 = 0
    n0::Float64 = nfunc(idx,raypath[1,1],raypath[1,2])
    x0::Float64 = raypath[1,1]*cos(raypath[1,2])
    y0::Float64 = raypath[1,1]*sin(raypath[1,2])
    
    for i = 2:size(raypath,1)
        x::Float64 = raypath[i,1]*cos(raypath[i,2])
        y::Float64 = raypath[i,1]*sin(raypath[i,2])
        ds::Float64 = hypot(x-x0,y-y0)
        pathlength += ds
        
        n::Float64 = nfunc(idx,raypath[i,1],raypath[i,2])
        action += ds*(n+n0)/2
        
        x0 = x
        y0 = y
        n0 = n
    end
    
    return (pathlength,action)
end

#Calculates the path integrals for specified segments
function getpathintegrals(raypath::Array{Float64,2},idx::RefractiveIndex,bounceindices::Array{Int64,1})
    #Initiate arrays
    bouncenum::Int64 = length(bounceindices)
    pathlengths::Array{Float64,1} = Array(Float64,bouncenum+1)
    actions::Array{Float64,1} = Array(Float64,bouncenum+1)
    
    #Start segment (before 1st bounce)
    pathlengths[1],actions[1] = getpathintegrals(raypath[1:bounceindices[1],:],idx)
    
    #Middle segment (between consecutive bounces)
    for i = 2:bouncenum
        pathlengths[i],actions[i] = 
            getpathintegrals(raypath[bounceindices[i-1]:bounceindices[i],:],idx)
    end
    
    #End segment (after last bounce)
    pathlengths[bouncenum+1],actions[bouncenum+1] =
        getpathintegrals(raypath[bounceindices[end]:end,:],idx)
    
    return pathlengths,actions
end