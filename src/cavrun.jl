#require("cavidxs.jl");
#require("cavbnds.jl");
require("cavsim.jl");
#require("cavcalc.jl");
require("cavplot.jl");
require("cavio.jl");

#The results are stored in a directory ../results with subdirectories for the cavity type (boundary shape and index distribution), each of which in turn has subdirectories for each set of cavity parameters.


#==============================================================================
#Helper functions


#Use PSS-spanning initial conditions
function getpssinitval(bnd::Boundary,idx::RefractiveIndex,symmetry::Int64=1)
    #Setup initial conditions
    thetares::Int64 = 25;#div(100,symmetry);
    sinchires::Int64 = 25;
    deviation::Int64 = 13; #coprime with and close to half of thetares
    theta0::Array{Float64,1} = linspace(0,2*pi/symmetry,thetares+1)[1:end-1];
    sinchi0::Array{Float64,1} = linspace(0,1,sinchires+3)[1:end-3];
    offset::Array{Float64,1} = linspace(0,sinchi0[deviation],thetares+1)[1:end-1];
    
    theta0val::Array{Float64,2} = ones(Float64,sinchires)*transpose(theta0);
    sinchi0val::Array{Float64,2} = mod(sinchi0*ones(Float64,1,thetares) + ones(Float64,sinchires)*transpose(offset),sinchi0[end]+sinchi0[2]) + sinchi0[2];
    
    r0val::Array{Float64,2},normangval::Array{Float64,2} = rsys(bnd,theta0val);
    phi0val::Array{Float64,2} = normangval + asin(sinchi0val);
    initval::Array{Float64,2} = [r0val[:] theta0val[:] phi0val[:]];
    
    return initval::Array{Float64,2};
end


#==============================================================================
#Run the simulator


#Analyze a single cavity and record results on disk
function runcavity(bnd::Boundary,idx::RefractiveIndex;symmetry::Int64=1,kwargs...)
    runcavity(bnd,idx,getpssinitval(bnd,idx,symmetry);kwargs...);
    return nothing;
end

function runcavity(bnd::Boundary,idx::RefractiveIndex,initval0::Array{Float64,2};maxuclst::Union(Float64,Array{Float64,1})=1.0,kwargs...)
    #Runs the ray tracer for a given cavity over multiple initial conditions, and
    #records the results of the simulation. Used mainly for obtaining the PSS and
    #calculating properties of the cavity.
    
    #maxuclst is the maximum unclusteredness score (as a scalar value or an array of 
    #values for each set of initial conditions) that a run will be recorded to disk, 
    #for better speed.
    assert((typeof(maxuclst) == Float64) || (length(maxuclst) == size(initval0,1)));
    
    #Check presence of results files
    cavdir::String = getcavdir(bnd,idx);
    if !isfile(joinpath(cavdir,"cavity.txt"))
        #Write header file with cavity information
        writehdr(cavdir,bnd,idx);
    else
        #Parse header file
        bnd1::Boundary,idx1::RefractiveIndex,cavlabel::String = parsehdr(cavdir);
        
        #Run only if identical cavity in header file description
        checkeq::Bool = ((names(bnd1)==names(bnd)) && (names(idx1)==names(idx)));
        if checkeq
            for field in names(bnd); checkeq &= (bnd1.(field)==bnd.(field)); end
            for field in names(idx); checkeq &= (idx1.(field)==idx.(field)); end
        end
        if !checkeq
            error("The current boundary and index labels are in use by a different cavity. This should not happen except in a hash collision.");
        end
    end
    
    #Save cavity plot image and get plot range
    cavplot = plt.figure(figsize=(5,5))[:gca](aspect="equal");
    plotbnd(bnd,axes=cavplot);
    plt.gcf()[:tight_layout]();
    range = plt.axis();
    if !isfile(joinpath(cavdir,"cavity.png"))
        plotidx(idx,bnd,axes=cavplot);
        plt.savefig(joinpath(cavdir,"cavity.png"));
    end
    plt.close();
    
    #Process init values
    if size(initval0,2) == 3
        initval::Array{Float64,2} = initval0;
    elseif size(initval0,2) == 2
        #Provided with [theta0 sin(chi0)]
        r0,normang = rsys(bnd,initval0[:,1]);
        phi0 = normang + asin(initval0[:,2]);
        initval = [r0 initval0[:,1] phi0];
    end
    
    #Begin runs for each set of initial conditions ----------------------------
    tic();
    runnum::Int64 = 0;
    for i = 1:size(initval,1)
        
        #Run ray tracer
        raypath::Array{Float64,2},bounces1::Array{Int32,1},bounces::Array{Float64,2} = rayevolve_gsl(bnd,idx,initval[i,:][:];kwargs...);
        
        #Calculate ray properties: pathlengths & actions between intervals, 
        #unclusteredness score, etc
        uclst::Float64, order::Int64 = unclusteredness(bounces);
        
        #Discard insufficiently clustered results
        #(Likely an uninteresting chaotic mode rather than one near a periodic mode)
        if (((typeof(maxuclst) == Float64) && (uclst > maxuclst)) 
            || ((typeof(maxuclst) == Array{Float64,1}) && (uclst > maxuclst[i])))
            continue;
        end
        
        #Continue processing results
        bouncenum::Int64 = size(bounces,1);
        pathlengths::Array{Float64,1} = Array(Float64,bouncenum+1);
        actions::Array{Float64,1} = Array(Float64,bouncenum+1);
        pathlengths[1], actions[1] = pathint(raypath[1:bounces1[1,1],:],idx);
        for j = 2:bouncenum
            pathlengths[j], actions[j] = pathint(raypath[bounces1[j-1,1]:bounces1[j,1],:],idx);
        end
        pathlengths[bouncenum+1], actions[bouncenum+1] = pathint(raypath[bounces1[end,1]:end,:],idx);
        
        #Find the next unassigned run number
        while isfile(joinpath(cavdir,"$(dec(runnum,4)).txt"))
            runnum += 1;
        end
        
        #Save ray plot
        rayplot = plt.figure(figsize=(5,5))[:gca](aspect="equal");
        plotrays(raypath,axes=rayplot);
        plt.axis(range);
        plt.gcf()[:tight_layout]();
        plt.axis("off");
        plt.savefig(joinpath(cavdir,"$(dec(runnum,4)).png"),transparent=true);
        plt.close();
        
        #Record results
        resultsfile::IOStream = open(joinpath(cavdir,"$(dec(runnum,4)).txt"),"a");
        writeresults(resultsfile,
            "init",initval[i,:][:], #line 2
            "pathlengths",pathlengths, #line 4
            "actions",actions, #line 6
            "bounces",bounces, #line 8
            "uclst order",(uclst,order)); #line 10
        
        if uclst < 0.15
            #Might meaningfully approximate a periodic mode
            
            #Calculate mode properties
            dwdm::Float64, omega0_s::Complex128, omega0_p::Complex128, phi_t::Array{Float64,1}, output_s::Array{Float64,1}, output_p::Array{Float64,1} = modeeval(order,pathlengths,actions,bounces,bnd,idx);
            
            #Record mode properties
            writeresults(resultsfile,
                "dwdm omega0_s omega0_p",(dwdm,omega0_s,omega0_p), #line 12
                "phi_t",phi_t, #line 14
                "output_s",output_s, #line 16
                "output_p",output_p); #line 18
        end
        close(resultsfile);
        
        #Display progress
        println("Completed $i of $(size(initval,1)) run(s).");
        
    end
    println("Cavity analyzed in $(toq()) seconds.");
    
    return nothing;
end


#Run cavity for mode position estimates without writing to disk
function getmodesest(bnd::Boundary,idx::RefractiveIndex;symmetry::Int64=1,kwargs...)
    return getmodesest(bnd,idx,getpssinitval(bnd,idx,symmetry);kwargs...);
end

function getmodesest(bnd::Boundary,idx::RefractiveIndex,initval0::Array{Float64,2};maxuclst::Union(Float64,Array{Float64,1})=1.0,maxrefines::Int64=5,tgtuclst::Float64=1e-4,kwargs...)
    
    assert((typeof(maxuclst) == Float64) || (length(maxuclst) == size(initval0,1)));
    
    #Initiate
    uclstval::Array{Float64,1} = Array(Float64,0);
    modebouncesval::Array{Array{Float64,2},1} = Array(Array{Float64,2},0);
    modeinitval::Array{Float64,2} = copy(initval0);
    
    count::Int64 = 0;
    unrefined::Array{Bool,1} = (uclstval .> tgtuclst);
    while (count <= maxrefines) && ((count == 0) || any(unrefined))
    
        #Process init values
        if size(modeinitval,2) == 3
            initval::Array{Float64,2} = modeinitval;
        elseif size(modeinitval,2) == 2
            #Provided with [theta0 sin(chi0)]
            initval = Array(Float64,size(modeinitval,1),3);
            initval[:,1],normang::Array{Float64,1} = rsys(bnd,modeinitval[:,1]);
            initval[:,2] = modeinitval[:,1];
            initval[:,3] = normang + asin(modeinitval[:,2]);
        end
        
        tic();
        for i = 1:size(initval,1)
            #Run ray tracer
            raypath::Array{Float64,2},bounces1::Array{Int32,1},bounces::Array{Float64,2} = rayevolve_gsl(bnd,idx,initval[i,:][:];kwargs...);
            
            #Calculate ray properties: pathlengths & actions between intervals, 
            #unclusteredness score, etc
            uclst::Float64, order::Int64 = unclusteredness(bounces);
            
            #Discard insufficiently clustered results
            #(Likely to be uninteresting mode rather than a periodic mode)
            if (((typeof(maxuclst) == Float64) && (uclst > maxuclst))
                || ((typeof(maxuclst) == Array{Float64,1}) && (uclst > maxuclst[i])))
                continue;
            end
            
            #Estimate mode bounce locations
            modenum::Int64 = 1;
            modebounces::Array{Float64,2} = getmodebounces(order,bounces);
            while modenum <= length(modebouncesval)
                if ((order == size(modebouncesval[modenum],1))
                    && (modedist(modebounces,modebouncesval[modenum]) < 0.002))
                    #Likely that both are approximating the same mode
                    if uclst < uclstval[modenum]
                        #Replace previous mode if better approximation
                        modebouncesval[modenum] = modebounces;
                        uclstval[modenum] = uclst;
                    end
                    break;
                end
                modenum += 1;
            end
            if modenum > length(modebouncesval)
                #New mode found
                push!(uclstval,uclst);
                push!(modebouncesval,modebounces);
            end
            
            #Display progress
            println("Completed $i of $(size(initval,1)) run(s).");
            
        end
        println("Cavity analyzed in $(toq()) seconds.");
        
        
        #Get the initial conditions for the mode estimates of the insufficiently 
        #clustered runs
        unrefined = (uclstval .> tgtuclst);
        modeinitval = Array(Float64,nnz(unrefined),2);
        for i = 1:nnz(unrefined)
            modeinitval[i,:] = modebouncesval[unrefined][i][1,:];
        end
        maxuclst = uclstval[unrefined];
        
        #Update mode collection
        uclstval = uclstval[!unrefined];
        modebouncesval = modebouncesval[!unrefined];
        
        count += 1;
    end
    
    return (uclstval,modebouncesval);
end


#Get cavity mode properties
function testcavity(bnd::Boundary,idx::RefractiveIndex,property::String="pss";kwargs...)
    #Refine mode estimates
    uclstval::Array{Float64,1}, mbval::Array{Array{Float64,2},1} = getmodesest(bnd,idx,maxuclst=0.6;kwargs...);
    
    #Get initial conditions
    minitval::Array{Float64,2} = Array(Float64,length(mbval),2);
    for i = 1:length(mbval)
        minitval[i,:] = mbval[i][1,:];
    end
    
    #Final run for data extraction
    if property == "pss"
        runcavity(bnd,idx,minitval,maxuclst=uclstval);
        showcavity(bnd,idx);
        return nothing;
        
    elseif property == "directionality"
        return NaN;
        
    elseif property == "mode stability"
        return NaN;
        
    else
        error("Property $property is not recognized.");
    end
end

