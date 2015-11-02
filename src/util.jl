# #############################################################################
# #############################################################################
#Helpful functions for computation


#Small nonnegative integer powers implemented in GSL
function gsl_pow_int(x::Float64,n::Int64)
    return ccall((:gsl_pow_int,"libgsl"),Cdouble,(Cdouble,Cint),x,convert(Int32,n))
end

#Squares and cubes are more common, simply implemented directly in Julia with shorter names, and also vectorized
sq{T<:Number}(x::T) = x*x::T
cb{T<:Number}(x::T) = x*x*x::T
@vectorize_1arg Number sq
@vectorize_1arg Number cb

#Custom Base64 representation of a number
const base64digits = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+-"
function dec2base64(n::Integer,pad::Int64=0)
    return lpad(base64digits[reverse(digits(n,64)+1)],pad,'0')
end

#Generic hashing of composite types based on typename, fields and values
function fieldhashing(obj)
    #Initiate hash with typename
    h::UInt64 = hash(summary(obj))
    for field in fieldnames(obj)
        #Use xor to combine hashes from all fields so that order does not matter
        #Hash the field name and value together as a tuple
        h $= hash((field,obj.(field)))
    end
    return h
end

#Generic equality of composite types based on typename, fields and values
function fieldequality(o1,o2)
    if typeof(o1) != typeof(o2)
        return false
    end
    for field in fieldnames(o1)
        if o1.(field) != o2.(field)
            return false
        end
    end
    return true
end

#Shrinking a matrix by deleting the last few rows or columns as much as possible without losing any nonzero value. Useful for matrices of coefficients, where position is important but zero values can be disregarded.
#This saves time interating through trailing zero coefficients, and also allows equivalent cavity objects with different coefficient matrix sizes to be hashed identically.
#Since this is run only at the creation of each instance of a Boundary or RefractiveIndex object, it does not affect computation time.
function coeffmatrix{T}(m::Array{T,1})
    I::Array{Int64,1},J::Array{Int64,1},V::Array{T,1} = findnz(cat(2,m))
    if isempty(I)
        return Array(T,0)
    else
        return full(sparsevec(I,V))[:]
    end
end
function coeffmatrix{T}(m::Array{T,2})
    I::Array{Int64,1},J::Array{Int64,1},V::Array{T,1} = findnz(m)
    if isempty(I)
        return Array(T,0,0)
    else
        return full(sparse(I,J,V))
    end
end

