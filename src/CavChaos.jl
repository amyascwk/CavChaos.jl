# #############################################################################
# #############################################################################
#Define module

module CavChaos
    
    #Pre-Julia v0.4.0 compatibility
    if VERSION < v"0.4-"
        UInt64 = Uint64
        AbstractString = String
        export UInt64, AbstractString
    else
        Union_type = Union
        Union(args...) = Union_type{args...}
        export Union
    end
    
    
    # #########################################################################
    #Cavity definition and miscellaneous layer
    
    #Internal functions used by other files
    include("util.jl")
    export dec2base64
    
    #Cavity boundary type definitions and functions
    include("boundary.jl")
    export Boundary, CircularBnd, EllipticalBnd, PolygonalBnd, FourierBnd
    export rfunc, rderiv, rderiv2, normang, rsys
    
    #Cavity refractive index distribution type definitions and functions
    include("refractive_index.jl")
    export RefractiveIndex, UniformIdx, GaussianIdx, TaylorIdx, ChebFourIdx, SumIdx
    export nfunc, nderiv
    
    #I/O functions for recording results
    include("io.jl")
    export hash
    export getresultsdir
    export writeresults, readresults
    export collectresultids
    
    
    # #########################################################################
    #Ray-solving layer
    
    #Wrapper function for ray evolution solver implemented in libcavchaos.so
    #C code in ray.c
    include("ray_solver.jl")
    export rayevolve_gsl
    
    #Coprocessing functions for analyzing ray data during simulation run,
    #(before recording results)
    include("ray_analysis.jl")
    export getpathintegrals
    
    #Ray plotting and image I/O functions
    include("ray_plot.jl")
    export plotbnd, plotidx, plotrays
    export getcavityimgfile, getrayimgfile
    export writecavityimg, writerayimg
    
    #Combined ray simulation
    include("ray_combined.jl")
    export run_rays, gen_pssinitarray, gen_initarray
    
    #Control file interface for single cavity run
    include("run_cavity.jl")
    export set_cavity_bnd!, get_cavity_bnd, set_cavity_idx!, get_cavity_idx
    export set_simulation_params!, get_simulation_param
    export set_init_params!, get_init_param
    export set_pssinitarray!, get_pssinitarray, set_initarray!, get_initarray
    export set_resultsroot!, get_resultsroot, get_resultsdir
    export set_results!, unset_results!, get_results
    export runcavity
    
    
    # #########################################################################
    #Mode-finding layer
    
    #Clustering algorithm
    include("clustering.jl")
    export findcluster
    
    #Wrapper function for mode finders implemented in libcavchaos.so
    #C code in mode.c
    include("mode_finder.jl")
    export findmode_gsl
    
    
    # #########################################################################
    #Optimization layer
    
    
end

