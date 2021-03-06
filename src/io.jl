# #############################################################################
# #############################################################################
#I/O methods

#This file contains functions for recording and reading the results of computations.

#   resultsdir = getresultsdir(simulation_params_hash::AbstractString,bnd::Boundary,idx::RefractiveIndex,resultsroot::AbstractString=".")
#   Returns the path <resultsdir> (as a subdirectory of <resultsroot>) to the uniquely assigned directory location for storing results of the simulation with parameters represented by the hash <simulation_params_hash>, for the cavity with boundary <bnd> and index distribution <idx>. Makes nonexisting directories via mkdir() whenever appropriate.

#   writeresults(resultsdir::AbstractString,resultid::Int64, label::AbstractString,data::Array)
#   writeresults(resultsdir::AbstractString,resultid::Int64, [label1::AbstractString,data1::Union{Integer,AbstractFloat},[label2::AbstractString,data2::Union{Integer,AbstractFloat},[...]]])
#   Writes the results corresponding to the <resultid>th computation of a specified simulation run, into files at the directory location <resultsdir> (as a subdirectory of <resultsroot>). For each label-data pair (<label>,<data>) in the variable argument list, the data are stored in appropriate formats (csv files) depending on their type, and systematically associated with their corresponding labels for later retrieval.
#   If multiple label-data argument pairs are present, either the data must all be of Integer or AbstractFloat types, so that they can be combined into a single file and identified by their position in the order. This is meant to reduce the clutter of multiple result files with only a single value stored in each. In other words, each call to writeresults should correspond to 1 result file (or directory for 3D and higher arrays).

#   data = readresults(resultsdir::AbstractString,resultid::Int64,label::AbstractString, T::Type=Float64)
#   Extracts the recorded data <data> associated to the label <label> for the <resultid>th computation of a specified simulation run, and stored in the directory <resultsdir>. If <data> is of type Number, attempts to reinterpret it as the type <T>.


# #############################################################################
# #############################################################################
#Initiate

#Dependencies
#require("util.jl")
#require("boundary.jl")
#require("refractive_index.jl")


# #############################################################################
# #############################################################################
#Results directory structure

#The result files are organized into directories first by cavity type (Boundary-RefractiveIndex constructor pairs), then by the specific parameters, then finally the actual control file that is executed to run the computations.

#Get cavity-specific directory
#Results directory is one more level deeper, separated by specific control file
function getcavitydir(bnd::Boundary,idx::RefractiveIndex,resultsroot::AbstractString=".")
    #Get cavity constructor names
    cavitytype::AbstractString = @sprintf("%s_%s",
                                          replace("$(typeof(bnd))",r"^CavChaos\.",""),
                                          replace("$(typeof(idx))",r"^CavChaos\.",""))
    #Get hash of specific cavity with parameters
    cavityhash = dec2base64(hash(bnd,hash(idx)))
    #Return path to cavity directory
    return joinpath(resultsroot,"results",cavitytype,"cav"*cavityhash)
end

#Get run-specific directory
function getresultsdir(run_params_hash::UInt64,bnd::Boundary,idx::RefractiveIndex,resultsroot::AbstractString=".";makedir::Bool=true)
    return joinpath(getcavitydir(bnd,idx,resultsroot),"run"*dec2base64(run_params_hash))
end


# #############################################################################
# #############################################################################
#Result file naming conventions

#Combine labels for multiple (nonarray) data in a single file
function combinelabels(labels...)
    return join(labels,"_")::AbstractString
end

#Split combined label into array of labels
function splitlabel(combinedlabel::AbstractString)
    return split(combinedlabel,'_')
end

#Get datafile filename from label (without extension)
#The resultid is included at the start of the filename if it is nonzero
function getresultfname(resultid::Int64,label::AbstractString)
    return (resultid == 0) ? label : @sprintf("%06d_%s",resultid,label)
end

#Extract the result id and labels from a datafile filename
function splitresultfname(fname::AbstractString)
    noextfname::AbstractString = splitext(fname)[1]
    #Check if datafile filename has a resultid at the start
    if ismatch(Regex("^\\d{6}_"),fname)
        return parse(noextfname[1:6]), splitlabel(noextfname[8:end])
    else
        return 0, splitlabel(noextfname)
    end
end

#Check if datafile filename has a given label
#Return (resultid,dataindex) pair, where dataindex is an indicator that is
#-1 if label was not found, 0 if label refers to data in whole file,
#and the actual index of the data otherwise
function checkfnameforlabel(fname::AbstractString,label::AbstractString)
    resultid::Int64,labels::Array{AbstractString,1} = splitresultfname(fname)
    index::Int64 = findfirst(labels,label)
    dataindex::Int64 = (index == 0) ? -1 : (length(labels) == 1) ? 0 : index
    return resultid, dataindex
end


# #############################################################################
# #############################################################################
#Result file I/O

#Write results for array data, with specified associated labels
#Data labels should not include '_' or '.', since these are used as separators of multiple data
function writeresults{T,N}(resultsdir::AbstractString,resultid::Int64,label::AbstractString,data::Array{T,N})
    if N <= 2
        #1D and 2D case
        #single csv file suffices for 1D or 2D arrays
        resultfile::AbstractString = joinpath(resultsdir,getresultfname(resultid,label)*".dat")
        writecsv(resultfile,data)
    else
        #3D and above case
        #need multiple files so store them in a directory
        resultdir::AbstractString = joinpath(resultsdir,getresultfname(resultid,label))
        mkdir(resultdir)
        datasizeN::Int64 = size(data,N)
        numberofdigits::Int64 = floor(log10(datasizeN))+1
        newdims = size(data)[1:end-1]
        for i=1:datasizeN
            #Recursively write results for each slice of data along the last dimension
            writeresults(resultdir,resultid,
                         #page labelling order will be opposite of dimension order
                         label*"."*dec(i,numberofdigits),
                         #the slice reshaped to (N-1)-dimensions
                         reshape(slicedim(data,N,i),newdims))
        end
    end
end

#Combine multiple single numerical data in a single file
function writeresults(resultsdir::AbstractString,resultid::Int64,args...)
    if !isempty(args)
        #Check arguments
        if isodd(length(args)); error("Argument list must contain label and data pairs."); end
        if !all(x -> typeof(x) <: Integer || typeof(x) <: AbstractFloat,
                args[2:2:end])
            error("Data types not suitable for recording in csv format. Only Integers and AbstractFloats accepted.")
        end
        #Write data
        writeresults(resultsdir,resultid,combinelabels(args[1:2:end]...),[args[2:2:end]...])
    end
end

#Read results
function readresults(resultsdir::AbstractString,resultid::Int64,label::AbstractString,T::Type=Float64)
    #First try simple data directory (for large dimensional arrays)
    resultdir::AbstractString = joinpath(resultsdir,getresultfname(resultid,label))
    if isdir(resultdir)
        #an N-dimensional array with N>2
        return readNdimcsv(resultdir,T)
    end
    
    #Otherwise try simple data file (for non-arrays or arrays of dimension < 3)
    if isfile(resultdir*".dat")
        #Matches exactly
        return readcsv(resultdir*".dat",T)
    end
    
    #Otherwise search for label in combined data files
    #Search for matches in list of files/directories
    for fname in readdir(resultsdir)
        resultid::Int64,dataindex::Int64 = checkfnameforlabel(fname,label)
        if dataindex > 0
            #found a file with multiple combined results,
            #one of which is the given label
            resultfile::AbstractString = joinpath(resultsdir,fname)
            if isfile(resultfile)
                entry = readcsv(resultfile)[dataindex]
                if typeof(entry) <: Number
                    return convert(T,entry)
                else
                    return entry
                end
            end
        end
    end
    
    #Otherwise throw an error
    error("The data associated to the label $label cannot be found for run ID $resultid in $resultsdir.")
end

#>> Read an N-dimensional array with N>2 from a directory
function readNdimcsv(resultpath::AbstractString,T::Type=Float64)
    if isfile(resultpath)
        return readcsv(resultpath,T)
    elseif isdir(resultpath)
        resultpartpaths::Array{AbstractString,1} = [joinpath(resultpath,fn) for fn in readdir(resultpath)]
        if isempty(resultpartpaths)
            #Return empty 2D array
            return Array(T,0,0)
        end
        #Get array of subarrays
        resultparts::Array = map(readNdimcsv,resultpartpaths,fill(T,length(resultpartpaths)))
        sizes = map(size,resultparts)
        #Check that all sizes match
        if all(sizes .== sizes[1])
            #Concatenate on next dimension
            return cat(ndims(resultparts[1])+1,resultparts...)
        else
            #size mismatch
            error("Size mismatch in subarrays when reading N-dimensional array from $resultpath.")
        end
    end
end


# #############################################################################
# #############################################################################
#Result collection

#Get list of resultids with a specific result recorded
function collectresultids(resultsdir::AbstractString,label::AbstractString)
    return map(first,
               filter(p -> p[2] != -1,
                      map(fn -> checkfnameforlabel(fn,label),
                          readdir(resultsdir))))
end
