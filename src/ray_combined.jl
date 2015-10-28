# #############################################################################
# #############################################################################
#Combined computation and analysis of rays in cavity

#This file contains the main function for running the ray solver on a cavity for multiple initial conditions, and storing the desired results on disk.

#   resultsdir = run_rays(bnd::Boundary,idx::RefractiveIndex;kwargs...)
#   This function runs the ray solvers and analysis on a given cavity with boundary <bnd> and index distribution <idx> on a number of initial conditions. The keyword arguments allow control over the following parameters regarding the initial conditions, the solver settings and results recording:
#
#>> Initial conditions
#   pssinitarray::Array{Float64,2} = Array(Float64,0,2)
#       array of initial conditions as (theta,chi) points on Poincare Surface of Section plot. Replaced with autogenerated points if empty.
#   initarray::Array{Float64,2} = Array(Float64,0,3)
#       array of initial conditions as (r,theta,phi) ray parameters. Replaced with equivalent values from <pssinitarray> if empty.
#   
#>> Initial conditions autogeneration parameters
#   symmetry::Int64 = 1
#       order of rotational symmetry of cavity to assume
#   thetares::Int64 = 25
#       number of theta values to sample
#   chires::Int64= 25
#       number of chi values to sample
#   
#>> Solver parameters
#   maxpathlength::Float64 = 200.0
#       maximum pathlength before halting solvers
#   maxbounces::Int64 = 500
#       maximum number of bounces before halting solvers
#   relativetolerance::Float64 = 1.0e-12
#       relative tolerance limits for solvers
#   absolutetolerance::Float64 = 1.0e-12
#       relative tolerance limits for solvers
#   
#>> Results root directory
#   resultsroot::String = "."
#       root directory to store result subdirectories
#   
#>> Results to record
#   record_pssinitarray::Bool = false
#   record_initarray::Bool = true
#   record_initialconditions::Bool = false
#   record_raypath::Bool = false
#   record_bounceindices::Bool = false
#   record_bouncepoints::Bool = true
#   record_cavityimage::Bool = true
#   record_rayimage::Bool = true
#   record_pathlengths::Bool = true
#   record_actions::Bool = true
#   record_modeproperties::Bool = false
#   record_farfield::Bool = false)
#       record value of parameter/result to disk


# #############################################################################
# #############################################################################
#Initiate

#Dependencies
#require("util.jl")
#require("boundary.jl")
#require("refractive_index.jl")
#require("solver.jl")
#require("io.jl")
#require("ray_plot.jl")
#require("ray_analysis.jl")


# #############################################################################
# #############################################################################
#Initial conditions generation

#Use PSS-spanning initial conditions
#Automatically generates an array of (theta,chi) pairs as rows, that covers as much of the Poincare Surface of Section plot up to a specified rotational symmetry.
function gen_pssinitarray(thetares::Int64=25,chires::Int64=25,symmetry::Int64=1)
    #Choose total number of points (see explanation below)
    interval::Int64 = int64(ceil(sqrt(thetares*chires)))
    N::Int64 = interval^2+3
    while gcd(interval,N) != 1
        N += 1
    end
    
    #Avoid repetition of values by generating N values of theta and N value of 
    #chi, then permute them relative to each other so that theta and chi values 
    #are matched to form points on a lattice spanning the PSS.
    #generate theta values spanning angle up to cavity rotational symmetry
    theta::Array{Float64,1} = linspace(0,2*pi/symmetry,N+1)[1:N]
    #generate chi values excluding 0 and pi/2
    chi::Array{Float64,1} = linspace(0,pi/2,N+2)[2:N+1]
    
    #In order to generate a lattice that is as evenly spread out as possible, can do 
    #this by stepping through appropriate intervals in mod <N>. By having N be in the 
    #form M^2 + P, where gcd(M,P) = 1, P < M, then lattice is fairly even with the P 
    #smallest values spread out as much as possible, in order to sample the low 
    #chi behavior at P points.
    return [theta chi[1+mod([(interval*(0:N-1))...],N)]]::Array{Float64,2}
end

#Converts a (theta,chi)-based set initial conditions to a (r,theta,phi)-based one
function pss2raycoord(bnd::Boundary,theta::Float64,chi::Float64)
    r::Float64,alpha::Float64 = rsys(bnd,theta)
    phi::Float64 = alpha + chi
    return (r,theta,phi)
end

#Converts a (r,theta,phi)-based initial conditions array from a given (theta,chi)-based array.
function pss2initarray(bnd::Boundary,pssinitarray::Array{Float64,2})
    count::Int64 = size(pssinitarray,1)
    initA::Array{Float64,2} = Array(Float64,count,3)
    for i=1:count
        initA[i,1],initA[i,2],initA[i,3] = pss2raycoord(bnd,pssinitarray[i,1],pssinitarray[i,2])
    end
    return initA
end

#Automatically generates an array of (r,theta,phi) as rows, that covers as much of the Poincare Surface of Section plot up to a specified rotational symmetry.
function gen_initarray(bnd::Boundary,thetares::Int64=25,chires::Int64=25,symmetry::Int64=1)
    pssinitarray::Array{Float64,2} = gen_pssinitarray(thetares,chires,symmetry)
    return pss2initarray(bnd,pssinitarray)
end


# #############################################################################
# #############################################################################
#Main ray simulation

#Run params hash
#Chunk all the simulation-relevant parameters together
#Then get a hash to use as a checksum
function get_runparamshash( initarray::Array{Float64,2},
                            maxpathlength::Float64,maxbounces::Int64,
                            relativetolerance::Float64,absolutetolerance::Float64,
                            record_pssinitarray::Bool,record_initarray::Bool,
                            record_initialconditions::Bool,record_raypath::Bool,
                            record_bounceindices::Bool,record_bouncepoints::Bool,
                            record_cavityimage::Bool,record_rayimage::Bool,
                            record_pathlengths::Bool,record_actions::Bool,
                            record_modeproperties::Bool,record_farfield::Bool)
    run_params = cell(3)
    run_params[1] = initarray
    run_params[2] = Number[maxpathlength,maxbounces,relativetolerance,absolutetolerance]
    run_params[3] = Bool[record_cavityimage,record_pssinitarray,record_initarray,record_initialconditions,record_raypath,record_bounceindices,record_bouncepoints,record_cavityimage,record_rayimage,record_pathlengths,record_actions,record_modeproperties,record_farfield]
    return hash(run_params)
end


#Actual run
function run_rays(  #Cavity parameters
                    bnd::Boundary,idx::RefractiveIndex;
                    
                    #Initial conditions
                    pssinitarray::Array{Float64,2} = Array(Float64,0,2),
                    initarray::Array{Float64,2} = Array(Float64,0,3),
                    
                    #Initial conditions generator parameters
                    symmetry::Int64 = 1,
                    thetares::Int64 = 25,
                    chires::Int64 = 25,
                    
                    #Solver parameters
                    maxpathlength::Float64 = 200.0,
                    maxbounces::Int64 = 500,
                    relativetolerance::Float64 = 1.0e-12,
                    absolutetolerance::Float64 = 1.0e-12,
                    
                    #Results root directory
                    resultsroot::String = ".",
                    
                    #Results to record
                    record_pssinitarray::Bool = false,
                    record_initarray::Bool = true,
                    record_initialconditions::Bool = false,
                    record_raypath::Bool = false,
                    record_bounceindices::Bool = false,
                    record_bouncepoints::Bool = true,
                    record_cavityimage::Bool = true,
                    record_rayimage::Bool = true,
                    record_pathlengths::Bool = true,
                    record_actions::Bool = true,
                    record_modeproperties::Bool = false,
                    record_farfield::Bool = false)
    
    tic()
    
    #Autogenerate initarray
    if isempty(initarray)
        if isempty(pssinitarray)
            pssinitarray = gen_pssinitarray(thetares,chires,symmetry)
        end
        initarray = pss2initarray(bnd,pssinitarray)
    end
    
    #Get results directory
    const run_params_hash::Uint64 = 
        get_runparamshash(  initarray,
                            maxpathlength,maxbounces,
                            relativetolerance,absolutetolerance,
                            record_pssinitarray,record_initarray,
                            record_initialconditions,record_raypath,
                            record_bounceindices,record_bouncepoints,
                            record_cavityimage,record_rayimage,
                            record_pathlengths,record_actions,
                            record_modeproperties,record_farfield)
    const resultsdir::String = getresultsdir(run_params_hash,bnd,idx,resultsroot)
    
    #Write cavity image and get plot range
    if record_cavityimage
        range::Array{Float64,1} = writecavityimg(resultsdir,bnd,idx)
    else
        #Set plot range to default
        cavityplot = prepareplot()
        plotbnd(bnd,axes=cavityplot)
        range = [plt.axis()...]
    end
    #Initialize ray plot line object to nothing
    rayline = nothing
    
    #Write initial conditions arrays
    if record_initarray
        writeresults(resultsdir,0,"initarray",initarray)
    end
    if record_pssinitarray
        writeresults(resultsdir,0,"pssinitarray",pssinitarray)
    end
    
    #Iterate through various initial conditions in parallel
    const num_of_runs::Int64 = size(initarray,1)
    @sync begin
        @parallel for resultid = 1:num_of_runs
            
            #Run ray tracer
            (raypath::Array{Float64,2},
             bounceindices::Array{Int64,1},
             bouncepoints::Array{Float64,2}) = 
                rayevolve_gsl(
                    bnd, idx,
                    initarray[resultid,:][:];
                    tmax = maxpathlength,
                    bouncemax = maxbounces,
                    reltol = relativetolerance,
                    abstol = absolutetolerance)
            
            #Analyze results
            if record_pathlengths || record_actions
                #Compute 
                pathlengths,actions = 
                    getpathintegrals(raypath,idx,bounceindices)
            end
            
            #Write results
            if record_initialconditions
                writeresults(resultsdir,resultid,
                             "initialconditions",initarray[resultid,:][:])
            end
            if record_raypath
                writeresults(resultsdir,resultid,"raypath",raypath)
            end
            if record_bounceindices
                writeresults(resultsdir,resultid,"bounceindices",bounceindices)
            end
            if record_rayimage
                if rayline == nothing
                    rayline = writerayimg(resultsdir,resultid,raypath,range)
                else
                    writerayimg(resultsdir,resultid,raypath,rayline)
                end
            end
            if record_bouncepoints
                writeresults(resultsdir,resultid,"bouncepoints",bouncepoints)
            end
            if record_pathlengths
                writeresults(resultsdir,resultid,"pathlengths",pathlengths)
            end
            if record_actions
                writeresults(resultsdir,resultid,"actions",actions)
            end
            
            #Mode detected! (not implemented yet)
            if false
                if record_modeproperties
                    
                end
                if record_farfield
                    
                end
            end
            
            #Display progress
            println("Completed #$resultid of $(num_of_runs) run(s).")
        end
    end
    
    #End
    plt.close()
    println("Cavity analyzed in $(toq()) seconds.")
    
    return resultsdir
end
