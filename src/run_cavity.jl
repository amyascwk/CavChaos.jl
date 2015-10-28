# #############################################################################
# #############################################################################
#Cavity analysis functions and global variables

#This file contains the main driving functions for evaluating a cavity, as well as the state variables and setting/getting functions for the simulation and cavity parameters. Typically, a control file will run on a multiprocess julia instance

#   shell-prompt> julia -p <num of parallel procs> <ctrlfile>

#The file starts with the declaration

#   @everywhere using CavChaos

#to initiate the module's methods into the namespace on every process. The following functions can then be used in a control file to set up a simulation run:

#>> Cavity properties

#   set_cavity_bnd!(bnd::Boundary)
#   Sets the desired cavity boundary to <bnd>.

#   bnd = get_cavity_bnd()
#   Returns the current cavity boundary.

#   set_cavity_idx!(idx::RefractiveIndex)
#   Sets the desired cavity index distribution to <idx>.

#   idx = get_cavity_idx()
#   Returns the current cavity index distribution.

#>> Solver parameters

#   set_solver_params!(;kwargs...)
#   Sets the solver parameter represented by keywords to the numerical value associated to it in <kwargs>. Appropriate keywords are: maxpathlength, maxbounces, relativetolerance, absolutetolerance.

#   value = get_solver_param(key::Union(AbstractString,Symbol))
#   Returns the current value of the solver parameter represented by the symbol or string <key>.

#>> Initial conditions

#   set_init_params!(;kwargs...)
#   Sets the initial conditions autogeneration parameter represented by keywords to the numerical value associated to it in <kwargs>. Appropriate keywords are: symmetry, thetares, chires.

#   value = get_init_param(key::Union(AbstractString,Symbol))
#   Returns the current value of the initial conditions autogeneration parameter represented by the symbol or string <key>.

#   set_pssinitarray!(A::Array{Float64,2})
#   Sets the custom initial conditions as an array <A> of initial (theta,chi) values as its rows. If nonempty, the autogeneration parameters will be ignored, and the initial conditions fed to the solver will be computed from <A> directly. This is useful for choosing points on the Poincare Surface of Section for performing simulation runs.

#   pssinitarray = get_pssinitarray()
#   Returns the current custom/autogenerated (theta,chi)-based initial conditions array <pssinitarray>.

#   set_initarray!(A::Array{Float64,2})
#   Sets the custom initial conditions as an array <A> of initial (r,theta,phi) values as its rows. If nonempty, the autogeneration parameters and the (theta,chi)-based custom initial conditions array will be ignored, and the initial conditions fed to the solver will be taken from <A> directly. This is useful for fine control over the initial ray position and direction when performing simulation runs.

#   initarray = get_initarray()
#   Returns the current custom/autogenerated (r,theta,phi)-based initial conditions array <initarray>.

#>> Results directories

#   set_resultsroot!(dir::AbstractString)
#   Sets the results root directory (in which a "results" subdirectory stores all output data) to <dir>.

#   get_resultsroot()
#   Returns the results root directory.

#   get_resultsdir()
#   Returns the configuration-specific results directory

#>> Results to record

#   set_results!(;kwargs...)
#   Sets the result represented by keywords to a boolean value associated to it in <kwargs>, indicating whether that result will be recorded. Appropriate keywords are: pssinitarray, initarray, initialconditions, raypath, bounceindices, bouncepoints, cavityimage, rayimage, pathlengths, actions, modeproperties, farfield.

#   set_results!(args...)
#   unset_results!(args...)
#   For a convenient syntax, if the above results keywords are passed as regular symbol or string arguments to set_results!, then they will be recorded. If passed to unset_results!, they will not.

#   boolean = get_results(key::Union(AbstractString,Symbol))
#   Returns whether the result represented by the keyword represented by the string or symbol <key> will be recorded.

#Finally, the actual run is started by a call to the function runcavity().


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
#require("ray_combined.jl")


# #############################################################################
# #############################################################################
#Settings as global state

#Cavity objects
cavity_bnd = FourierBnd()
cavity_idx = UniformIdx()

#Solver parameters
solver = Dict{Symbol,Number}()
solver[:maxpathlength] = 200.0
solver[:maxbounces] = 500
solver[:relativetolerance] = 1.0e-12
solver[:absolutetolerance] = 1.0e-12

#Initial conditions array
init = Dict{Symbol,Number}()
init[:symmetry] = 1
init[:thetares] = 25
init[:chires]= 25
pssinitarray = Array(Float64,0,2)
initarray = Array(Float64,0,3)

#Root directory to store results
resultsroot = "."

#Results to record
results = Dict{Symbol,Bool}()
results[:pssinitarray] = false
results[:initarray] = true
results[:initialconditions] = false
results[:raypath] = false
results[:bounceindices] = false
results[:bouncepoints] = true
results[:cavityimage] = true
results[:rayimage] = true
results[:pathlengths] = true
results[:actions] = true
results[:modeproperties] = false
results[:farfield] = false


# #############################################################################
# #############################################################################
#Setters and unsetters and getters

#Cavity parameters
function set_cavity_bnd!(bnd::Boundary)
    global cavity_bnd = bnd
    return nothing
end

function get_cavity_bnd()
    return cavity_bnd::Boundary
end

function set_cavity_idx!(idx::RefractiveIndex)
    global cavity_idx = idx
    return nothing
end

function get_cavity_idx()
    return cavity_idx::RefractiveIndex
end

#Solver parameters
function set_solver_params!(;kwargs...)
    for i=1:length(kwargs)
        if haskey(solver,kwargs[i][1])
            solver[kwargs[i][1]] = kwargs[i][2]
        else
            error("No solver parameter $(kwargs[i][1]).")
        end
    end
    return nothing
end

function get_solver_param(key::Union(AbstractString,Symbol))
    return solver[symbol(key)]::Number
end

#Initial conditions parameters
function set_init_params!(;kwargs...)
    for i=1:length(kwargs)
        if haskey(init,kwargs[i][1])
            init[kwargs[i][1]] = kwargs[i][2]
        else
            error("No initial conditions parameter $(kwargs[i][1]).")
        end
    end
    return nothing
end

function get_init_param(key::Union(AbstractString,Symbol))
    return init[symbol(key)]::Number
end

function set_pssinitarray!(A::Array{Float64,2})
    if size(A,2) != 2
        error("PSS initial conditions array must have (theta,chi) pairs as rows.")
    end
    global pssinitarray = A
    return nothing
end

function get_pssinitarray()
    return pssinitarray::Array{Float64,2}
end

function set_initarray!(A::Array{Float64,2})
    if size(A,2) != 3
        error("Initial conditions array must have (r,theta,phi) values as rows.")
    end
    global initarray = A
    return nothing
end

function get_initarray()
    return initarray::Array{Float64,2}
end

#Root directory to store results
function set_resultsroot!(directory::AbstractString)
    global resultsroot = directory
    return nothing
end

function get_resultsroot()
    return resultsroot::AbstractString
end

#Results directory for the current configuration
function get_resultsdir()
    #Autogenerate initarray
    initarray = get_initarray()
    if isempty(initarray)
        pssinitarray = get_pssinitarray()
        if isempty(pssinitarray)
            pssinitarray = gen_pssinitarray(get_init_param(:thetares),
                                            get_init_param(:chires),
                                            get_init_param(:symmetry))
        end
        initarray = pss2initarray(get_cavity_bnd(),pssinitarray)
    end
    
    #Get run parameter hash
    run_params_hash::UInt64 = 
        get_runparamshash(  initarray,
                            get_solver_param(:maxpathlength),
                            get_solver_param(:maxbounces),
                            get_solver_param(:relativetolerance),
                            get_solver_param(:absolutetolerance),
                            get_results(:pssinitarray),get_results(:initarray),
                            get_results(:initialconditions),get_results(:raypath),
                            get_results(:bounceindices),get_results(:bouncepoints),
                            get_results(:cavityimage),get_results(:rayimage),
                            get_results(:pathlengths),get_results(:actions),
                            get_results(:modeproperties),get_results(:farfield))
    
    return getresultsdir(   run_params_hash,get_cavity_bnd(),get_cavity_idx(),
                            get_resultsroot(),makedir=false)
end

#Results parameters
function set_results!(;kwargs...)
    for i=1:length(kwargs)
        if haskey(results,kwargs[i][1])
            results[kwargs[i][1]] = kwargs[i][2]
        else
            error("No result named $(kwargs[i][1]).")
        end
    end
    return nothing
end

function set_results!(args...)
    set_results!(;[(symbol(arg),true) for arg in args]...)
end

function unset_results!(args...)
    set_results!(;[(symbol(arg),false) for arg in args]...)
end

function get_results(key::Union(AbstractString,Symbol))
    return results[symbol(key)]::Bool
end


# #############################################################################
# #############################################################################
#Main driver

#Profile a single cavity by running multiple rays in it
function runcavity()
    
    run_rays(   #Cavity parameters
                get_cavity_bnd(),get_cavity_idx();
                
                #Initial conditions
                pssinitarray = get_pssinitarray(),
                initarray = get_initarray(),
                
                #Initial conditions generator parameters
                symmetry = get_init_param(:symmetry),
                thetares = get_init_param(:thetares),
                chires= get_init_param(:chires),
                
                #Solver parameters
                maxpathlength = get_solver_param(:maxpathlength),
                maxbounces = get_solver_param(:maxbounces),
                relativetolerance = get_solver_param(:relativetolerance),
                absolutetolerance = get_solver_param(:absolutetolerance),
                
                #Results root directory
                resultsroot = get_resultsroot(),
                
                #Results to record
                record_pssinitarray = get_results(:pssinitarray),
                record_initarray = get_results(:initarray),
                record_initialconditions = get_results(:initialconditions),
                record_raypath = get_results(:raypath),
                record_bounceindices = get_results(:bounceindices),
                record_bouncepoints = get_results(:bouncepoints),
                record_cavityimage = get_results(:cavityimage),
                record_rayimage = get_results(:rayimage),
                record_pathlengths = get_results(:pathlengths),
                record_actions = get_results(:actions),
                record_modeproperties = get_results(:modeproperties),
                record_farfield = get_results(:farfield),
            )
    
end

