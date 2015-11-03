# #############################################################################
# #############################################################################
#Define module

module CavChaos
    
    #Pre-Julia v0.4.0 compatibility
    if VERSION < v"0.4-"
        #Getting PyPlot functions to work using plt[:func] syntax
        import Base.getindex
        getindex(m::Module,s::Symbol) = m.(s)

        #New type names
        typealias UInt64 Uint64
        typealias AbstractFloat FloatingPoint
        typealias AbstractString String

        #Rounding function replacements
        import Base.ceil
        ceil(T::DataType,x::Integer) = convert(T,ceil(x))
        ceil(T::DataType,x::Integer,y::Integer) = convert(T,ceil(x,y))
        ceil(T::DataType,args...) = convert(T,ceil(args...))
        function cld{T1<:Real,T2<:Real}(x::T1,y::T2)
            d,v = divrem(x,y)
            return (v <= 0) ? d : d + 1
        end
        
        #Introspection
        fieldnames = names

        #Hashing
        import Base.hash
        hash(h1::UInt64,h2::UInt64) = hash((h1,h2))
        
        export UInt64, AbstractString, fieldnames
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
    export getresultsdir
    export writeresults, readresults
    export collectresultids
    
    
    # #########################################################################
    #Ray-solving layer
    
    #Wrapper function for ray evolution solver implemented in libcavchaos.so
    #C code in ray.c
    include("ray_solver.jl")
    export rayevolve_gsl
    
    #Clustering algorithm (specific to periodic points on PSS)
    include("clustering.jl")
    export findcluster
    
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
    export set_solver_params!, get_solver_param
    export set_init_params!, get_init_param
    export set_pssinitarray!, get_pssinitarray, set_initarray!, get_initarray
    export set_resultsroot!, get_resultsroot, get_resultsdir
    export set_results_params!, unset_results_params!, get_results_param
    export run_cavity
    
    
    # #########################################################################
    #Mode-finding layer
    
    #PSS plotting and image I/O functions
    include("mode_plot.jl")
    export plotpss
    export getpssimgfile, writepssimg
    
    #Wrapper function for mode finders implemented in libcavchaos.so
    #C code in mode.c
    include("mode_finder.jl")
    export findmode_gsl
    
    
    # #########################################################################
    #Optimization layer
    
    
end

