#require("cavidxs.jl");
#require("cavbnds.jl");
require("cavsim.jl");
reload("cavcalc.jl");
require("cavplot.jl");

#The results are stored in a directory ../results with subdirectories for the cavity type (boundary shape and index distribution), each of which in turn has subdirectories for each set of cavity parameters.

#==============================================================================
#Parse files


#Generate Julia-readable output
jr(x) = string(x);
function jr{T}(x::Array{T,1})
    s::String = replace("$T[$(string(x)[1:end-1])]","\n",",");
    return s;
end
function jr{T}(x::Array{T,2})
    if string(x) == ""; return "Array($T,$(size(x)))"; end
    s::String = replace("$T[$(string(x)[1:end-1])]","\n",";");
    s = replace(s,"\t"," ");
    return s;
end
function jr(x::String)
    return "\"$x\"";
end


#Get cavity information from header file
function parseheader(filepath::String)
    hdrfile::IOStream = open(filepath,"r");
    filelines::Array{String,1} = readlines(hdrfile);
    close(hdrfile);
    
    cavlabel::String = "";
    bnd::Boundary = CircularBnd();
    idx::RefractiveIndex = UniformIdx();
    linenum::Int64 = findfirst(filelines,"[Boundary]\n")+1;
    count::Int64 = 0;
    while (linenum <= length(filelines)) && (filelines[linenum][1] != '[')
        if filelines[linenum][1] != '#'
            eqnpos::Int64 = findfirst(filelines[linenum],'=')
            if eqnpos != 0
                #Read boundary field values
                field::Symbol = parse(filelines[linenum][1:eqnpos-1]);
                if field == :boundary
                    bnd = eval(parse(filelines[linenum][eqnpos+1:end-1]));
                else
                    bnd.(field) = eval(parse(filelines[linenum][eqnpos+1:end-1]));
                end
            end
            if (count >= 2) && (filelines[linenum][1] != '\n')
                str1 = replace(filelines[linenum],"_","\\_");
                cavlabel = "$(cavlabel)$(str1)";
            end
        end
        linenum += 1;
        count += 1;
    end
    
    linenum = findfirst(filelines,"[Index]\n")+1;
    count = 0;
    while (linenum <= length(filelines)) && (filelines[linenum][1] != '[')
        if filelines[linenum][1] != '#'
            eqnpos = findfirst(filelines[linenum],'=')
            if eqnpos != 0
                #Read boundary field values
                field = parse(filelines[linenum][1:eqnpos-1]);
                if field == :index
                    idx = eval(parse(filelines[linenum][eqnpos+1:end-1]));
                else
                    idx.(field) = eval(parse(filelines[linenum][eqnpos+1:end-1]));
                end
            end
            if (count >= 2) && (filelines[linenum][1] != '\n')
                str1 = replace(filelines[linenum],"_","\\_");
                cavlabel = "$(cavlabel)$(str1)";
            end
        end
        linenum += 1;
        count += 1;
    end
    return (bnd,idx,cavlabel);
    
end


#==============================================================================
#Run the simulator


#Analyze a single cavity
function runcavity(bnd::Boundary,idx::RefractiveIndex,initval::Array{Float64,2};kwargs...)
    #Runs the ray tracer for a given cavity over multiple initial conditions, and
    #records the results of the simulation. Used mainly for obtaining the PSS and
    #calculating properties of the cavity.
    
    #Get directory for cavity type and parameters
    cavdir::String = joinpath("..","results");
    if !isdir(cavdir)
        mkdir(cavdir);
    end
    #Cavity type
    cavtype::String = "$(summary(bnd)[1:end-3])_$(summary(idx)[1:end-3])";
    cavdir = joinpath(cavdir,"$cavtype");
    if !isdir(cavdir)
        mkdir(cavdir);
    end
    #Cavity parameters (specified in cavity boundary and index names)
    runname::String = "$(bnd.name)_$(idx.name)";
    cavdir = joinpath(cavdir,"$runname");
    if !isdir(cavdir)
        mkdir(cavdir);
    end
    
    #Check presence of important files
    if !isfile(joinpath(cavdir,"cavity.txt"))
        #Write header file with cavity information
        hdrfile::IOStream = open(joinpath(cavdir,"cavity.txt"),"w");
        println(hdrfile,"# PLEASE DO NOT EDIT ANY FILE IN THIS DIRECTORY");
        println(hdrfile,"# Cavity boundary ################################");
        println(hdrfile,"[Boundary]");
        println(hdrfile,"boundary = $(summary(bnd))()");
        for field in names(bnd)[1:end-3]
            println(hdrfile,"$(field) = $(jr(bnd.(field)))");
        end
        println(hdrfile,"\n# Cavity index distribution ######################")
        println(hdrfile,"[Index]");
        println(hdrfile,"index = $(summary(idx))()");
        for field in names(idx)[1:end-3]
            println(hdrfile,"$(field) = $(jr(idx.(field)))");
        end
        close(hdrfile);
    else
        #Parse header file
        bnd1::Boundary,idx1::RefractiveIndex,cavlabel::String = parseheader(joinpath(cavdir,"cavity.txt"));
        
        #Run only if identical cavity in header file description
        checkeq::Bool = ((names(bnd1)==names(bnd)) && (names(idx1)==names(idx)));
        if checkeq
            for field in names(bnd)[1:end-3]
                checkeq &= (bnd1.(field) == bnd.(field));
            end
            for field in names(idx)[1:end-3]
                checkeq &= (idx1.(field) == idx.(field));
            end
        end
        if !checkeq
            error("The current boundary and index labels are in use by a different cavity.");
        end
    end
    
    #Save cavity plot image and get plot range
    cavplot = plt.figure(figsize=(5,5))[:gca](aspect="equal");
    plotbnd(bnd,axes=cavplot);
    plotidx(idx,bnd,axes=cavplot);
    plt.gcf()[:tight_layout]();
    range = plt.axis();
    if !isfile(joinpath(cavdir,"cavity.png"))
        plt.savefig(joinpath(cavdir,"cavity.png"));
    end
    plt.close();
    
    #Run and Record
    tic();
    runnum::Int64 = 0;
    for i = 1:size(initval,1)
        #Find the next unassigned run number
        while isfile(joinpath(cavdir,"$(dec(runnum,4)).txt"))
            runnum += 1;
        end
        
        #Run ray tracer
        raypath::Array{Float64,2},bounces::Array{Float64,2} = rayevolve_gsl(initval[i,1],initval[i,2],initval[i,3],bnd=bnd,idx=idx,tmax=200.0,kwargs...);
        
        #Save ray plot
        rayplot = plt.figure(figsize=(5,5))[:gca](aspect="equal");
        plotrays(raypath,axes=rayplot);
        plt.axis(range);
        plt.axis("off");
        plt.gcf()[:tight_layout]();
        plt.savefig(joinpath(cavdir,"$(dec(runnum,4)).png"),transparent=true);
        plt.close();
        
        #Calculate pathlengths and actions
        bouncenum::Int64 = size(bounces,1);
        pathlengths::Array{Float64,1} = Array(Float64,bouncenum+1);
        actions::Array{Float64,1} = Array(Float64,bouncenum+1);
        pathlengths[1], actions[1] = pathint(raypath[1:int64(bounces[1,1]),:],idx);
        for j = 2:bouncenum
            pathlengths[j], actions[j] = pathint(raypath[int64(bounces[j-1,1]):int64(bounces[j,1]),:],idx);
        end
        pathlengths[bouncenum+1], actions[bouncenum+1] = pathint(raypath[int64(bounces[end,1]):end,:],idx);
        bounces[:,1] = raypath[int64(bounces[:,1]),2]; #replace index with theta value
        
        #Record results
        resultsfile::IOStream = open(joinpath(cavdir,"$(dec(runnum,4)).txt"),"a");
        println(resultsfile,"# init");
        println(resultsfile,"$(jr(initval[i,:][:]))");
        println(resultsfile,"# pathlengths");
        println(resultsfile,"$(jr(pathlengths))");
        println(resultsfile,"# actions");
        println(resultsfile,"$(jr(actions))");
        println(resultsfile,"# bounces");
        println(resultsfile,"$(jr(bounces))");
        
        #Calculate mode properties
        uclst::Float64, order::Int64 = interact[:unclusteredness](bounces);
        if uclst < 0.15
            #Close to a periodic mode
            dwdm::Float64, omega0_s::Complex128, omega0_p::Complex128, phi_t::Array{Float64,1}, output_s::Array{Float64,1}, output_p::Array{Float64,1} = modeeval(order,pathlengths,actions,bounces,bnd,idx);
            
            #Record mode properties
            println(resultsfile,"# uclst order");
            println(resultsfile,"($uclst, $order)");
            println(resultsfile,"# dwdm omega0_s omega0_p");
            println(resultsfile,"($dwdm, $omega0_s, $omega0_p)");
            println(resultsfile,"# phi_t");
            println(resultsfile,"$(jr(phi_t))");
            println(resultsfile,"# output_s");
            println(resultsfile,"$(jr(output_s))");
            println(resultsfile,"# output_p");
            println(resultsfile,"$(jr(output_p))");
        end
        
        close(resultsfile);
        
        #Display progress
        println("Completed $i of $(size(initval,1)) run(s).");
    end
    println("Cavity analyzed in $(toq()) seconds.");
    
    return nothing;
end
#Use PSS-spanning initial conditions
function runcavity(bnd::Boundary,idx::RefractiveIndex;symmetry::Int64=1,kwargs...)
    #Setup initial conditions
    thetares = 25;#div(100,symmetry);
    sinchires = 27;
    theta0::Array{Float64,1} = linspace(0,2*pi/symmetry,thetares);
    sinchi0::Array{Float64,1} = linspace(0,1,sinchires)[2:end-1];
    offset = linspace(0,sinchi0[1],length(theta0)+1)[1:end-1];
    
    theta0val = ones(Float64,length(sinchi0))*transpose(theta0);
    sinchi0val = mod(sinchi0*ones(Float64,1,length(theta0)) + ones(Float64,length(sinchi0))*transpose(offset),1.0);
    
    r0val::Array{Float64,2},normangval::Array{Float64,2} = bnd.sys(theta0val);
    phi0val = normangval + asin(sinchi0val);
    initval = [r0val[:] theta0val[:] phi0val[:]];
    runcavity(bnd,idx,initval,kwargs...);
    return nothing;
    
end


#==============================================================================
#Display results


#Present PSS
function showcavity(cavdir::String)
    #Present the results of computation on a cavity
    
    #Read cavity information
    if (!isfile(joinpath(cavdir,"cavity.txt")))
        error("No ray tracing data found.");
    end
    bnd::Boundary,idx::RefractiveIndex,cavlabel::String = parseheader(joinpath(cavdir,"cavity.txt"));
    
    #Read simulation results
    init::Array{Float64,1} = Array(Float64,3);
    bouncesval::Array{Array{Float64,2},1} = Array(Array{Float64,2},0);
    
    #Read simulation results
    runnum::Int64 = 0;
    while isfile(joinpath(cavdir,"$(dec(runnum,4)).txt"))
        resultsfile::IOStream = open(joinpath(cavdir,"$(dec(runnum,4)).txt"),"r");
        filelines::Array{String,1} = readlines(resultsfile);
        close(resultsfile);
        push!(bouncesval,eval(parse(filelines[8])));
        if runnum == 0 
            init = eval(parse(filelines[2]));
        end
        runnum += 1;
    end
    
    #Generate Plots -------------------------------------------------------
    #Farfield emissions plot
    fig2 = plt.figure("farfield");
    
    fig = plt.figure("PSS explorer");
    fig[:subplots_adjust](left=0.05,bottom=0.1,right=0.95,top=0.9,wspace=0.05,hspace=0.3)
    cavinfo = plt.subplot2grid((5,2),(0,0),title="\\Large{\\textbf{Cavity: $(summary(bnd)[1:end-3])\\_$(summary(idx)[1:end-3])/$(bnd.name)\\_$(idx.name)}}");
    rayplot = plt.subplot2grid((5,2),(1,0),rowspan=4,aspect="equal");
    pssplot = plt.subplot2grid((5,2),(0,1),rowspan=5,title="Poincare Surface of Section");
    fig[:canvas][:mpl_connect]("pick_event",interact["onpick"]);
    #fig[:canvas][:mpl_connect]("button_press_event",interact["onpress"]);
    
    #Ray plot (Initialize with data for first run)
    initlabel = "Run 0000: \$(r_0,\\theta_0,\\phi_0)\$ = ($(init[1]), $(init[2]), $(init[3]))\n\$(\\theta_0,\\sin\\chi_0)\$ = ($(bouncesval[1][1]), $(sin(bouncesval[1][2])))";
    rayplot[:set](title="Ray Plot",xlabel=initlabel,xticks=[],yticks=[]);
    cavimg = plt.imread(joinpath(cavdir,"cavity.png")); #Background image
    rayplot[:imshow](cavimg);
    rayimg = plt.imread(joinpath(cavdir,"0000.png")); #Ray foreground image
    rayplot[:imshow](rayimg,axes=rayplot);
    
    #Cavity label
    if cavlabel[end] == '\n'; cavlabel = cavlabel[1:end-1]; end
    modelabel::String = "Click any set of points on the PSS plot\nto the right.";
    cavinfo[:text](0,1,cavlabel,ha="left",va="top",size=12);
    cavinfo[:text](0.5,1,modelabel,ha="left",va="top",size=12);
    mpl[:rc]("text",usetex=true);
    cavinfo[:axis]("off");
    
    #PSS plot
    plotpss(bouncesval,bnd,idx,axes=pssplot);
    
    plt.show();
    # ---------------------------------------------------------------------
    
    return nothing;
end

