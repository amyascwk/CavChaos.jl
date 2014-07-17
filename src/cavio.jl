#require("cavbnds.jl");
require("cavidxs.jl");


#==============================================================================
#Results files I/O


#Cavity hash function
function cavhash(bnd::Boundary,idx::RefractiveIndex)
    fieldhashes::Array{Uint64,1} = Array(Uint64,0);
    for field in names(bnd)
        push!(fieldhashes,hash(bnd.(field)));
    end
    for field in names(idx)
        push!(fieldhashes,hash(idx.(field)));
    end
    return hex(hash(fieldhashes))::String;
end


#Get directory for cavity type and parameters
function getcavdir(bnd::Boundary,idx::RefractiveIndex)
    cavdir::String = "..";
    cavtype::String = "$(summary(bnd)[1:end-3])_$(summary(idx)[1:end-3])";
    for subdir::String in ["results",cavtype,cavhash(bnd,idx)]
        cavdir = joinpath(cavdir,subdir);
        if !isdir(cavdir); mkdir(cavdir); end
    end
    return cavdir::String;
end


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


#Write header (cavity information) file
function writehdr(cavdir::String,bnd::Boundary,idx::RefractiveIndex)
    hdrfile::IOStream = open(joinpath(cavdir,"cavity.txt"),"w");
    println(hdrfile,"# Cavity boundary ################################");
    println(hdrfile,"[Boundary]");
    println(hdrfile,"$(summary(bnd))()");
    for field in names(bnd)
        println(hdrfile,"$(field) = $(jr(bnd.(field)))");
    end
    println(hdrfile,"\n# Cavity index distribution ######################")
    println(hdrfile,"[Index]");
    println(hdrfile,"$(summary(idx))()");
    for field in names(idx)
        println(hdrfile,"$(field) = $(jr(idx.(field)))");
    end
    close(hdrfile);
end


#Get cavity information from header file
function parsehdr(cavdir::String)
    hdrfile = open(joinpath(cavdir,"cavity.txt"),"r");
    filelines::Array{String,1} = readlines(hdrfile);
    cavlabel::String = "";
    
    linenum::Int64 = findfirst(filelines,"[Boundary]\n")+1;
    bnd::Boundary = eval(parse(filelines[linenum][1:end-1]));
    while (linenum <= length(filelines)) && (filelines[linenum][1] != '[')
        if filelines[linenum][1] != '#'
            eqnpos::Int64 = findfirst(filelines[linenum],'=')
            if eqnpos != 0
                #Read boundary field values
                field::Symbol = parse(filelines[linenum][1:eqnpos-1]);
                bnd.(field) = eval(parse(filelines[linenum][eqnpos+1:end-1]));
                line::String = replace(filelines[linenum],"_","\\_");
                if length(line) > 50; line = "$(line[1:46]) ...\n"; end
                cavlabel = "$(cavlabel)$(line)";
            end
        end
        linenum += 1;
    end
    
    linenum = findfirst(filelines,"[Index]\n")+1;
    idx::RefractiveIndex = eval(parse(filelines[linenum][1:end-1]));
    while (linenum <= length(filelines)) && (filelines[linenum][1] != '[')
        if filelines[linenum][1] != '#'
            eqnpos = findfirst(filelines[linenum],'=')
            if eqnpos != 0
                #Read boundary field values
                field = parse(filelines[linenum][1:eqnpos-1]);
                idx.(field) = eval(parse(filelines[linenum][eqnpos+1:end-1]));
                line = replace(filelines[linenum],"_","\\_");
                if length(line) > 50; line = "$(line[1:46]) ...\n"; end
                cavlabel = "$(cavlabel)$(line)";
            end
        end
        linenum += 1;
    end
    close(hdrfile);
    
    return (bnd,idx,cavlabel);
end


#Write results
function writeresults(resultsfile::IOStream,label::String,data,args...)
    println(resultsfile,"# $label");
    println(resultsfile,jr(data));
    if length(args) > 0
        writeresults(resultsfile,args...);
    end
    return nothing;
end