# #############################################################################
# #############################################################################
#PSS clustering methods

#This file contains functions to quantify and analyze the degree of clustering of points on a Poincare Surface of Section, which would indicate the presence of a nearby periodic mode.

#Finding the order of periodicity of a series of bounces from their positions on the PSS is a constrained clustering problem, where the data is already ordered such that the clusters are sampled periodically, so points from the same cluster are separated by uniform intervals equal to the number of clusters.


# #############################################################################
# #############################################################################
#Distance definitions on the Poincare Surface of Section

#Metric on the PSS
function pssmetric(dtheta::Float64,dchi::Float64)
    return hypot(2*dchi,sin(dtheta/2))
end

#Distance function between 2 points on the PSS
#Column vectors
function pssdist(p1::Array{Float64,1},p2::Array{Float64,1})
    assert(length(p1) == length(p2) == 2)
    return pssmetric(p1[1]-p2[1],p1[2]-p2[2])
end
#Row vectors
function pssdist(p1::Array{Float64,2},p2::Array{Float64,2})
    assert(size(p1,1) == size(p2,1) == 1)
    return pssdist(p1[:],p2[:])
end


# #############################################################################
# #############################################################################
#Algorithm to estimate number of clusters

#This algorithm relies on the following properties of clusters in the specific context of k-periodic bounce points:

#   - points an interval of k apart in the bounce sequence are separated by a very 
#     small PSS distance.
#   - due to the fractal nature of the PSS, the clusters of bounce points themselves 
#     appear in clusters of a divisor of k.

#This means that one can progressively test increasing intervals until a sharp drop in PSS distance indicates a periodic cluster, and subsequently increment the test intervals by the number of clusters, since k should be some multiple of it. The details of the algorithm are as follows:

#   1) Distance between pairs of points a fixed interval apart in the bounce sequence 
#      is estimated by the mean of such pairwise distances from all such pairs 
#      available in the bounce data. This approximates the distance between their 
#      corresponding clusters, so a sudden low value indicates that each pair belongs #      to the same cluster.
#   2) The cluster count is initially estimated to be 1
#   3) The cluster distance for various intervals of bounce points is computed, with 
#      the interval incremented by the cluster count, since that should be a divisor 
#      of the "actual" order of periodicity.
#   4) If a sharp drop is detected below a threshhold proportion of the 
#      typical cluster distances, then it is likely that there is a clustering 
#      present. The cluster count is now estimated to be the current test interval.
#   5) Repeat (3) and (4) until the test intervals is too high relative to the total 
#      number of bounce points for meaningful clusters to form. If the cluster count 
#      remains 1, then it is inconclusive whether there is any clustering present at 
#      all.

#Mean PSS distance of pairs of bouncepoints a specified interval apart
function periodic_pssdist(bncpts::Array{Float64,2},order::Int64,samplesize::Int64=0)
    const maxsamplesize::Int64 = size(bncpts,1)-order
    if samplesize == 0 || samplesize > maxsamplesize
        samplesize = maxsamplesize
    end
    return mean([pssdist(bncpts[i,:],bncpts[i+order,:]) for i=1:samplesize])
end

#Estimate the number of clusters of bouncepoints
function findcluster(bncpts::Array{Float64,2};threshhold::Float64=0.1,samplesize::Int64=30)
    #Initiate
    clustcount::Int64 = 1
    testinterval::Int64 = 2
    clustsepsum::Float64 = periodic_pssdist(bncpts,1,samplesize)
    clustsepcount::Int64 = 1
    
    const maxtestinterval::Int64 = size(bncpts,1)-samplesize
    while(testinterval < maxtestinterval)
        #Get mean PSS distance for test interval
        testdist::Float64 = periodic_pssdist(bncpts,testinterval,samplesize)
        
        #Test against threshhold
        if testdist < threshhold*clustsepsum/clustsepcount
            #Sharp drop detected
            #cluster count likely to be multiple of test interval
            clustcount = testinterval
            #Reset cluster distance data
            clustsepsum = testdist
            clustsepcount = 1
        else
            #No significant drop detected
            #Collect cluster separation data
            clustsepsum += testdist
            clustsepcount += 1
        end
        
        #Increment test interval by the current cluster count estimate
        testinterval += clustcount
    end
    
    #Return best estimate of number of clusters
    return clustcount, clustsepsum/clustsepcount
end

