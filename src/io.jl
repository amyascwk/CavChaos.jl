# #############################################################################
# #############################################################################
#I/O methods

#This file contains functions for recording and reading the results of computations.

#   resultsdir = getresultsdir(simulation_params_hash::AbstractString,bnd::Boundary,idx::RefractiveIndex,resultsroot::AbstractString=".")
#   Returns the path <resultsdir> (as a subdirectory of <resultsroot>) to the uniquely assigned directory location for storing results of the simulation with parameters represented by the hash <simulation_params_hash>, for the cavity with boundary <bnd> and index distribution <idx>. Makes nonexisting directories via mkdir() whenever appropriate.

#   writeresults(resultsdir::AbstractString,resultid::Int64, [label1::AbstractString,data1,[label2::AbstractString,data2,[...]]])
#   Writes the results corresponding to the <resultid>th computation of a specified simulation run, into files at the directory location <resultsdir> (as a subdirectory of <resultsroot>). For each label-data pair (<label>,<data>) in the variable argument list, the data are stored in appropriate formats (csv files) depending on their type, and systematically associated with their corresponding labels for later retrieval.
#   If multiple label-data argument pairs are present, either the data must all be of array-type, which are then each stored in individual files, or the data must all NOT be of array-type, and they are combined into a single file and identified by their position in the order. This is to reduce the clutter of multiple result files with only a single value stored in each.

#   data = readresults(resultsdir::AbstractString,resultid::Int64,label::AbstractString, T::Type=Float64)
#   Extracts the recorded data <data> associated to the label <label> for the <resultid>th computation of a specified simulation run, and stored in the directory <resultsdir>. If <data> is of type Number, attempts to reinterpret it as the type <T>.


# #############################################################################
# #############################################################################
#Initiate

#Dependencies
#require("util.jl")
#require("boundary.jl")
#require("refractive_index.jl")
import Base.hash


# #############################################################################
# #############################################################################
#Hashing cavity objects

#Extend base hash function with method to hash equivalent Boundary and RefractiveIndex objects identically.
function hash(cavobj::Union(Boundary,RefractiveIndex))
    infoarray = Any[summary(cavobj)]
    for field in names(cavobj)
        push!(infoarray,cavobj.(field))
    end
    return hash(infoarray)
end


# #############################################################################
# #############################################################################
#Results directory structure

#The result files are organized into directories first by cavity type (Boundary-RefractiveIndex constructor pairs), then by the specific parameters, then finally the actual control file that is executed to run the computations.

#Get cavity-specific directory
#Results directory is one more level deeper, separated by specific control file
function getcavitydir(bnd::Boundary,idx::RefractiveIndex,resultsroot::AbstractString=".")
    #Initiate cavity directory as module root directory
    cavitydir::AbstractString = resultsroot
    #Get cavity constructor names
    cavitytype::AbstractString = @sprintf("%s_%s",summary(bnd),summary(idx))
    #Get hash of specific cavity with parameters
    cavityhash = dec2base64(hash((bnd,idx)))
    #Make path to cavity directory
    for subdir::AbstractString in ["results",cavitytype,"cav"*cavityhash]
        cavitydir = joinpath(cavitydir,subdir)
        if !isdir(cavitydir); mkdir(cavitydir); end
    end
    return cavitydir
end

#Get run-specific directory
function getresultsdir(run_params_hash::UInt64,bnd::Boundary,idx::RefractiveIndex,resultsroot::AbstractString=".";makedir::Bool=true)
    #Get directory name
    resultsdir::AbstractString = joinpath(getcavitydir(bnd,idx,resultsroot),"run"*dec2base64(run_params_hash))
    #Make directory
    if makedir && !isdir(resultsdir); mkdir(resultsdir); end
    return resultsdir
end


# #############################################################################
# #############################################################################
#Result file naming conventions

#Combine labels for multiple (nonarray) data in a single file
function combinelabels(label::AbstractString,labels...)
    return label*"_"*combinelabels(labels...)::AbstractString
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
        return int64(noextfname[1:6]), splitlabel(noextfname[8:end])
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

#Write results for different types of data, with specified associated labels
#Data labels should not include '_' or '.', since these are used as separators of multiple data

#>> Do nothing in empty case
writeresults(resultsdir::AbstractString,resultid::Int64) = nothing

#>> Write arrays directly
#   (1D and 2D case)
function writeresults{T}(resultsdir::AbstractString,resultid::Int64,label::AbstractString,data::Union(Array{T,1},Array{T,2}),args...)
    #single csv file suffices for 1D or 2D arrays
    resultfile::AbstractString = joinpath(resultsdir,getresultfname(resultid,label)*".dat")
    writecsv(resultfile,data)
    #Continue with other label-data argument pairs
    writeresults(resultsdir,resultid,args...)
end

#   (N dimensional case)
function writeresults{T,N}(resultsdir::AbstractString,resultid::Int64,label::AbstractString,data::Array{T,N},args...)
    #this is run when N>2
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
    #Continue with other label-data argument pairs
    writeresults(resultsdir,resultid,args...)
end

#>> Combine non-arrays into a single file
#   (assumes that all data arguments are non-arrays)
function writeresults(resultsdir::AbstractString,resultid::Int64,args...)
    data_array::Array{Any,1} = [args[2:2:end]...]
    for item in data_array
        assert(!(typeof(item) <: Array))
    end
    writeresults(resultsdir,resultid,combinelabels(args[1:2:end]...),data_array)
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
