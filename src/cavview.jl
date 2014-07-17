#require("cavbnds.jl");
#require("cavidxs.jl");
#require("cavsim.jl");
#require("cavcalc.jl");
require("cavplot.jl");
require("cavio.jl");


#==============================================================================
#Display results


#Plot Poincare Surface of Section (Memory efficient version)
function plotpss(cavdir::String;axes=plt.gca())
    #Setup plot
    plt.plot([-pi,pi],[0,0],"k",axes=axes);
    plt.axis([-pi,pi,-0.6,1]);
    plt.xlabel("\$\\displaystyle\\theta\$", fontsize=18);
    plt.ylabel("\$\\displaystyle\\sin\\chi\$", fontsize=18);
    
    #Parse results files
    runnum::Int64 = 0;
    while isfile(joinpath(cavdir,"$(dec(runnum,4)).txt"))
        resultsfile::IOStream = open(joinpath(cavdir,"$(dec(runnum,4)).txt"),"r");
        filelines::Array{String,1} = readlines(resultsfile);
        close(resultsfile);
        
        uclst::Float64 = eval(parse(filelines[10]))[1];
        bounces::Array{Float64,2} = eval(parse(filelines[8]));
        if runnum == 0 
            init = eval(parse(filelines[2]));
            plt.plot(bounces[:,1],sin(bounces[:,2]),".",ms=8,mec="k",mfc="r",mew=1.5,picker=true,zorder=1000,label="(\"0000\",$uclst)",axes=axes)
        else
            color = interact[:colorfunc](uclst);
            plt.plot(bounces[:,1],sin(bounces[:,2]),".",color=color,markersize=2,picker=true,label="(\"$(dec(runnum,4))\",$uclst)",axes=axes);
        end
        
        if mod(runnum,100) == 0; plt.draw(); end
        runnum += 1;
    end
    
    return nothing;
end


#Present PSS
function showcavity(bnd::Boundary,idx::RefractiveIndex)
    #Automatically find cavity directory
    cavdir::String = getcavdir(bnd,idx);
    if !isdir(cavdir); error("Cavity directory $cavdir not found."); end
    showcavity(cavdir);
    return nothing;
end

function showcavity(cavdir::String)
    #Present the results of computation on a cavity
    
    #Read cavity information
    if !isfile(joinpath(cavdir,"cavity.txt"))
        error("No ray tracing data found.");
    end
    bnd::Boundary,idx::RefractiveIndex,cavlabel::String = parsehdr(cavdir);
    
    #Read first result
    resultsfile::IOStream = open(joinpath(cavdir,"0000.txt"),"r");
    filelines::Array{String,1} = readlines(resultsfile);
    close(resultsfile);
    init = eval(parse(filelines[2]));
    bounces::Array{Float64,2} = eval(parse(filelines[8]));
    
    #Generate Plots -------------------------------------------------------
    plt.ion();
    mpl[:rc]("text",usetex=true);
    
    #Farfield Emissions figure
    fig2 = plt.figure("Far Field Emissions (still testing)");
    ffplot = plt.gca();
    
    #PSS Explorer figure
    fig = plt.figure("PSS Explorer",figsize=(16,9));
    fig[:subplots_adjust](left=0.05,bottom=0.1,right=0.95,top=0.9,wspace=0.05,hspace=0.3)
    cavtitle::String = replace(cavdir[12:end],"_","\\_");
    cavinfo = plt.subplot2grid((5,2),(0,0),title="\\Large{\\textbf{Cavity: $(cavtitle)}}");
    rayplot = plt.subplot2grid((5,2),(1,0),rowspan=4,aspect="equal");
    pssplot = plt.subplot2grid((5,2),(0,1),rowspan=5,title="Poincare Surface of Section (loading...)");
    fig[:canvas][:mpl_connect]("pick_event",interact["onpick"]);
    fig[:canvas][:mpl_connect]("button_press_event",interact["onpress"]);
    fig[:canvas][:mpl_connect]("scroll_event",interact["onscroll"]);
    
    #Ray plot (Initialize with data for first run)
    initlabel = "Run 0000: (\$r_0\$, \$\\theta_0\$, \$\\phi_0\$) = ($(init[1]), $(init[2]), $(init[3]))\n(\$\\theta_0\$, \$\\sin\\chi_0\$) = ($(bounces[1]), $(sin(bounces[2])))";
    rayplot[:set](title="Ray Plot",xlabel=initlabel,xticks=[],yticks=[]);
    cavimg = plt.imread(joinpath(cavdir,"cavity.png")); #Background image
    rayplot[:imshow](cavimg);
    rayimg = plt.imread(joinpath(cavdir,"0000.png")); #Ray foreground image
    rayplot[:imshow](rayimg,axes=rayplot);
    
    #Cavity label
    if cavlabel[end] == '\n'; cavlabel = cavlabel[1:end-1]; end #Strip trailing newline
    modelabel::String = "Click any set of points on the PSS plot\nto the right.";
    cavinfo[:text](0,0.5,cavlabel,ha="left",va="center",size=12);
    cavinfo[:text](0.5,0.5,modelabel,ha="left",va="center",size=12);
    cavinfo[:axis]("off");
    
    #PSS plot
    plotpss(bnd,idx,cavdir,axes=pssplot);
    pssplot[:findobj](mpl[:patches]["Rectangle"])[1][:set_label]("10");
    plt.title("Poincare Surface of Section");
    
    plt.show();
    plt.ioff();
    # ---------------------------------------------------------------------
    
    return nothing;
end


#Interactively list existing results
function listcavities()
    cavdir::String = joinpath("..","results");
    tmppath::String, tmpfile::IOStream = mktemp();
    println(tmpfile,"   | Cavity name     | Boundary  | Index     | Parameters");
    println(tmpfile,repeat("-",100));
    
    #Compile list of all cavity results available
    paths::Array{String,1} = Array(String,0);
    count::Int64 = 1;
    strm1,proc = readsfrom(`ls $cavdir`);
    while process_running(proc); run(`sleep 0.01`); end
    for cavtype::String in readlines(strm1)
    
        cavtype = cavtype[1:end-1];
        if !isdir(joinpath(cavdir,cavtype)); continue; end
        strm2,proc = readsfrom(`ls $(joinpath(cavdir,cavtype))`);
        while process_running(proc); run(`sleep 0.01`); end
        for cavity::String in readlines(strm2)
        
            cavity = cavity[1:end-1];
            if !isdir(joinpath(cavdir,cavtype,cavity)); continue; end
            if !isfile(joinpath(cavdir,cavtype,cavity,"cavity.txt")); continue; end
            
            bnd,idx,cavlabel = parsehdr(joinpath(cavdir,cavtype,cavity));
            
            push!(paths,joinpath(cavtype,cavity));
            print(tmpfile,"$(lpad(string(count),3,' ')): $(rpad(cavity,16,' '))  $(rpad(summary(bnd)[1:end-3],10,' '))  $(rpad(summary(idx)[1:end-3],10,' '))");
            values::String = "";
            for field::Symbol in names(bnd)
                value::String = jr(bnd.(field));
                num = findfirst(value,'[');
                if num != 0; value = value[num:end]; end
                if length(value) > 18; value = "$(value[1:14]) ..."; end
                values = "$(values)  $(value)";
            end
            for field in names(idx)
                value = jr(idx.(field));
                num = findfirst(value,'[');
                if num != 0; value = value[num:end]; end
                if length(value) > 18; value = "$(value[1:14]) ..."; end
                values = "$(values)  $(value)";
            end
            if length(values) > 70; values = values[1:70]; end
            print(tmpfile,"$(values)\n");
            count += 1;
            
        end
        close(strm2);
        println(tmpfile,"-");
    end
    close(strm1);
    close(tmpfile);
    
    #Show list
    run(`less $tmppath`);
    
    #Run input loop
    while true
        print("Enter task (\"l\" to list cavities again, \"h\" for help, \"q\" to quit): ");
        mode = readline(STDIN);
        if mode == "\n"
            #nothing
        elseif mode[1] == 'q'
            break;
            
        elseif mode[1] == 'h'
            #help
            println("<num>\t\tView information about cavity <num>.");
            println("pss <num>\tShow PSS plot of recorded data for cavity <num>.");
            println("show <num>\tShow image of cavity <num>.");
            println("rm <num>\tPermanently delete all data for cavity <num>.\n");
        
        elseif mode[1] == 'l'
            #list cavities again
            listcavities();
            break;
            
        elseif ismatch(r"^\s*\d+\s*$",mode[1:end-1]) && (1 <= int64(mode[1:end-1]) < count)
            #View cavity information
            run(`less $(joinpath(cavdir,paths[int64(mode[1:end-1])],"cavity.txt"))`);
            
        elseif ((mode[1:3] == "pss") && ismatch(r"^\s*\d+\s*$",mode[4:end-1])
                && (1 <= int64(mode[4:end-1]) < count))
            #Show PSS plot
            showcavity(joinpath(cavdir,paths[int64(mode[5:end-1])]));
            
        elseif ((mode[1:4] == "show") && ismatch(r"^\s*\d+\s*$",mode[5:end-1])
                && (1 <= int64(mode[5:end-1]) < count)
                && isfile(joinpath(cavdir,paths[int64(mode[5:end-1])],"cavity.png")))
            #Show cavity
            imgdata = plt.imread(joinpath(cavdir,paths[int64(mode[5:end-1])],"cavity.png"));
            plt.imshow(imgdata);
            plt.axis("off");
            plt.show();
            
        elseif ((mode[1:2] == "rm") && ismatch(r"^\s*\d+\s*$",mode[3:end-1])
                && (1 <= int64(mode[3:end-1]) < count))
            #Remove cavity results
            run(`/bin/rm -Ir $(joinpath(cavdir,paths[int64(mode[3:end-1])]))`);
            
        end
    end
    return nothing;
end


