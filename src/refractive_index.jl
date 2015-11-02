# #############################################################################
# #############################################################################
#Cavity refractive index distribution

#This file contains functions for specifying the refractive index distribution of the cavity, by subtypes of the abstract type RefractiveIndex. The following are the available constructors, which can take both positional or keyword arguments:

#   UniformIdx(base::Float64=2.0)
#   constructs a uniform distribution of the specified base value (default 2.0).

#   GaussianIdx(base::Float64=2.0,center::Array{Float64,1}=[0.0,0.0], amplitude::Float64=-0.2,spread::Array{Float64,1}=[0.5,0.5],rotation::Float64=0.0)
#   constructs a distribution with a Gaussian perturbation of specified center position (default at (x,y)=(0.0,0.0)), amplitude (default -0.2), standard deviation of spread in orthogonal directions (default dx=dy=0.5), and angle of rotation (default=0.0), over a specified base value (default 2.0).

#   TaylorIdx(coeff::Array{Float64,2}=[2.0 0.0;0.0 0.4])
#   constructs a general distribution specified by the Taylor expansion coefficients with respect to the x and y coordinates (default is n = 2.0 + 0.4*x*y).

#   ChebFourIdx(coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array(Float64,(0,0)))
#   constructs a general differentiable distribution specified by the coefficients for the cross terms of the odd Chebyshev expansion of the radial dependence and the Fourier expansion of the angular dependence. The even terms of the Chebyshev expansion are omitted because they are nonzero at the origin, which would make the distribution discontinuous at the origin if there is also a theta dependence. 
#   (Default is n = 2.0*T_0 + 0.1*T_1(r)*cos(theta) - 0.1*T_1(r)*cos(2*theta) - 0.1*T_3(r)*cos(theta) + 0.1*T_3(r)*cos(2*theta))

#   SumIdx(idxs...)
#   constructs a distribution that is the sum of the distributions specified by the arguments. Useful for multiple independent Gaussian perturbations.

#Information about the index distribution in polar coordinates can be extracted by the following functions, which dispatch according to the distribution function subtype to the appropriate calculations:

#   n = nfunc(idx::RefractiveIndex,r,theta)
#   returns the value of the refractive index distribution <idx> at the polar position (<r>,<theta>). <r> can be a Float64 or Array{Float64,N}, but <theta> must follow its type and dimensions, and the result <n> will also follow its type and dimensions.

#   n,dr_n,dtheta_n = nderiv(idx::RefractiveIndex,r,theta)
#   returns a tuple (<n>,<dr_n>,<dtheta_n>) of the value <n> and partial derivatives <dr_n> and <dtheta_n> of the refractive index distribution <idx> at the polar position (<r>,<theta>). <r> can be a Float64 or Array{Float64,N}, but <theta> must follow its type and dimensions, and each entry in the tuple result will also follow its type and dimensions. There are no functions to compute the derivatives alone since the calculations involved almost always requires the refractive index value to also be computed.

#   nderiv!(idx::RefractiveIndex,r::Float64,theta::Float64,result::Array{Float64,1})
#   This is the in-place version of nderiv storing the results as the first three entries of the array <result>.


# #############################################################################
# #############################################################################
#Initialize

#Dependencies
#require("util.jl")
#require("boundary.jl")
import Base: ==, hash

#General index distribution type
abstract RefractiveIndex


# #############################################################################
# #############################################################################
#Hashing and equality

#Extend base hash function with method to hash equivalent RefractiveIndex objects identically.
hash(idx::RefractiveIndex) = fieldhashing(idx)

#Equality is equality of type and field values
==(idx1::RefractiveIndex,idx2::RefractiveIndex) = fieldequality(idx1,idx2)


# #############################################################################
# #############################################################################
#Functions for generic index distribution

#These functions are vectorized/in-place versions of nfunc and nderiv for a general RefractiveIndex. The implementations for specific index distribution subtypes, if present, are slightly more efficient.

#Index value ------------------------------------------------------------------
function nfunc{N}(idx::RefractiveIndex,r::Array{Float64,N},theta::Array{Float64,N})
    #assert(size(r)==size(theta))
    #Initialize array
    n::Array{Float64,N} = similar(r)
    #Set value to nfunc output
    for i=1:length(r)
        n[i] = nfunc(idx,r[i],theta[i])
    end
    return n::Array{Float64,N}
end

#Index value and derivatives --------------------------------------------------
function nderiv{N}(idx::RefractiveIndex,r::Array{Float64,N},theta::Array{Float64,N})
    #assert(size(r)==size(theta))
    #Initialize arrays 
    n::Array{Float64,N} = similar(r)
    dr_n::Array{Float64,N} = similar(r)
    dtheta_n::Array{Float64,N} = similar(r)
    #Set values to nderiv output
    for i=1:length(theta)
        n[i],dr_n[i],dtheta_n[i] = nderiv(idx,r[i],theta[i])
    end
    return (n,dr_n,dtheta_n)
end

#In-place storing version of nderiv
function nderiv!(idx::RefractiveIndex,r::Float64,theta::Float64,result::Array{Float64,1})
    result[1],result[2],result[3] = nderiv(idx,r,theta)
    return nothing
end


# #############################################################################
# #############################################################################
#List of primitive refractive index distributions

# #############################################################################
#Uniform

#Uniform distribution (base value still needs to be specified to determine losses at boundary)
type UniformIdx <: RefractiveIndex
    #Base value
    base::Float64
    
    #Construct by fieldname
    UniformIdx(;base::Float64=2.0) = new(base)
    #Construct by position
    UniformIdx(base::Float64=2.0) = new(base)
end
        
#Index value ------------------------------------------------------------------
nfunc(idx::UniformIdx,r::Float64,theta::Float64) = idx.base::Float64
function nfunc{N}(idx::UniformIdx,r::Array{Float64,N},theta::Array{Float64,N})
    #assert(size(r)==size(theta))
    n::Array{Float64,N} = fill(idx.base,size(r))
    return n::Array{Float64,N}
end

#Index value and derivatives --------------------------------------------------
nderiv(idx::UniformIdx,r::Float64,theta::Float64) = (idx.base,0.0,0.0)
function nderiv{N}(idx::UniformIdx,r::Array{Float64,N},theta::Array{Float64,N})
    #assert(size(r)==size(theta))
    return (fill(idx.base,size(r)), zeros(Float64,size(r)), zeros(Float64,size(r)))
end


# #############################################################################
#Gaussian

#A Gaussian perturbation to a uniform base value.
type GaussianIdx <: RefractiveIndex
    #Base value
    base::Float64
    #Center position of the Gaussian pertubation (Array of x and y coordinate)
    center::Array{Float64,1}
    #Amplitude of the perturbation
    amplitude::Float64
    #The spreading of the Gaussian perturbation (Array of x and y stddev values)
    spread::Array{Float64,1}
    #Angle of rotation
    rotation::Float64
    
    #Construct by fieldname
    function GaussianIdx(;base::Float64=2.0,center::Array{Float64,1}=[0.0,0.0],amplitude::Float64=-0.2,spread::Array{Float64,1}=[0.5,0.5],rotation::Float64=0.0)
        return new(base,center,amplitude,spread,rotation)
    end
    #Construct by position
    function GaussianIdx(base::Float64=2.0,center::Array{Float64,1}=[0.0,0.0],amplitude::Float64=-0.2,spread::Array{Float64,1}=[0.5,0.5],rotation::Float64=0.0)
        return new(base,center,amplitude,spread,rotation)
    end
end
    
#Index value ------------------------------------------------------------------
function nfunc(idx::GaussianIdx,r::Float64,theta::Float64)
    x::Float64 = r*cos(theta)
    y::Float64 = r*sin(theta)
    dxr::Float64 = (x-idx.center[1])*cos(idx.rotation) - (y-idx.center[2])*sin(idx.rotation)
    dyr::Float64 = (y-idx.center[2])*cos(idx.rotation) + (x-idx.center[1])*sin(idx.rotation)
    n::Float64 = idx.base + idx.amplitude*exp(-0.5*(sq(dxr/idx.spread[1]) + sq(dyr/idx.spread[2])))
    return n
end

#Index value and derivatives --------------------------------------------------
function nderiv(idx::GaussianIdx,r::Float64,theta::Float64)
    x::Float64 = r*cos(theta)
    y::Float64 = r*sin(theta)
    dxr::Float64 = (x-idx.center[1])*cos(idx.rotation) - (y-idx.center[2])*sin(idx.rotation)
    dyr::Float64 = (y-idx.center[2])*cos(idx.rotation) + (x-idx.center[1])*sin(idx.rotation)
    gaussian::Float64 = idx.amplitude*exp(-0.5*(sq(dxr/idx.spread[1]) + sq(dyr/idx.spread[2])))
    
    n::Float64 = idx.base + gaussian
    dr_n::Float64 = -(cos(theta+idx.rotation)*dxr/sq(idx.spread[1]) + sin(theta+idx.rotation)*dyr/sq(idx.spread[2]))*gaussian
    dtheta_n::Float64 = ((y*cos(idx.rotation)+x*sin(idx.rotation))*dxr/sq(idx.spread[1]) + (-x*cos(idx.rotation)+y*sin(idx.rotation))*dyr/sq(idx.spread[2]))*gaussian
    return (n,dr_n,dtheta_n)
end


# #############################################################################
#Taylor

#Generic distribution generated by Tayor expansion along the x and y directions about the origin
type TaylorIdx <: RefractiveIndex
    #coeff is a 2D array where the (i,j)th element is the coefficient of the x^(i-1)*y^(j-1) term in the Taylor expansion about the origin of the index distribution. coeff needs only to be large enough to include the significant coefficients
    coeff::Array{Float64,2}
    
    #Construct by fieldname
    TaylorIdx(;coeff::Array{Float64,2}=[2.0 0.0;0.0 0.4]) = new(coeffmatrix(coeff))
    #Construct by position
    TaylorIdx(coeff::Array{Float64,2}=[2.0 0.0;0.0 0.4]) = new(coeffmatrix(coeff))
end

#Index value ------------------------------------------------------------------
function nfunc(idx::TaylorIdx,r::Float64,theta::Float64)
    x::Float64 = r*cos(theta)
    y::Float64 = r*sin(theta)
    n::Float64 = 0.0
    for i = 1:size(idx.coeff,1)
        for j = 1:size(idx.coeff,2)
            n += idx.coeff[i,j]*gsl_pow_int(y,i-1)*gsl_pow_int(x,j-1)
        end
    end
    return n::Float64
end

function nfunc{N}(idx::TaylorIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta))
    n::Array{Float64,N} = similar(r)
    for k = 1:length(r)
        x::Float64 = r[k]*cos(theta[k])
        y::Float64 = r[k]*sin(theta[k])
        n[k] = 0.0
        for i = 1:size(idx.coeff,1)
            for j = 1:size(idx.coeff,2)
                n[k] += idx.coeff[i,j]*gsl_pow_int(y,i-1)*gsl_pow_int(x,j-1)
            end
        end
    end
    return n::Array{Float64,N}
end

#Index value and derivatives --------------------------------------------------
function nderiv(idx::TaylorIdx,r::Float64,theta::Float64)
    x::Float64 = r*cos(theta)
    y::Float64 = r*sin(theta)
    n::Float64 = 0.0
    dr_n::Float64 = 0.0
    dtheta_n::Float64 = 0.0
    for i = 1:size(idx.coeff,1)
        for j = 1:size(idx.coeff,2)
            n += idx.coeff[i,j]*gsl_pow_int(y,i-1)*gsl_pow_int(x,j-1)
            if i >= 2
                dr_n += idx.coeff[i,j]*(i-1)*gsl_pow_int(y,i-2)*gsl_pow_int(x,j-1)*sin(theta)
                dtheta_n += idx.coeff[i,j]*(i-1)*gsl_pow_int(y,i-2)*gsl_pow_int(x,j)
            end
            if j >= 2
                dr_n += idx.coeff[i,j]*gsl_pow_int(y,i-1)*(j-1)*gsl_pow_int(x,j-2)*cos(theta)
                dtheta_n -= idx.coeff[i,j]*gsl_pow_int(y,i)*(j-1)*gsl_pow_int(x,j-2)
            end
        end
    end
    return (n,dr_n,dtheta_n)
end

function nderiv{N}(idx::TaylorIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta))
    n::Array{Float64,N} = similar(r)
    dr_n::Array{Float64,N} = similar(r)
    dtheta_n::Array{Float64,N} = similar(r)
    for k = 1:length(r)
        x::Float64 = r[k]*cos(theta[k])
        y::Float64 = r[k]*sin(theta[k])
        n[k] = 0.0
        dr_n[k] = 0.0
        dtheta_n[k] = 0.0
        for i = 1:size(idx.coeff,1)
            for j = 1:size(idx.coeff,2)
                n[k] += idx.coeff[i,j]*gsl_pow_int(y,i-1)*gsl_pow_int(x,j-1)
                if i >= 2
                    dr_n[k] += idx.coeff[i,j]*(i-1)*gsl_pow_int(y,i-2)*gsl_pow_int(x,j-1)*sin(theta[k])
                    dtheta_n[k] += idx.coeff[i,j]*(i-1)*gsl_pow_int(y,i-2)*gsl_pow_int(x,j)
                end
                if j >= 2
                    dr_n[k] += idx.coeff[i,j]*gsl_pow_int(y,i-1)*(j-1)*gsl_pow_int(x,j-2)*cos(theta[k])
                    dtheta_n[k] -= idx.coeff[i,j]*gsl_pow_int(y,i)*(j-1)*gsl_pow_int(x,j-2)
                end
            end
        end
    end
    return (n,dr_n,dtheta_n)
end


# #############################################################################
#Chebyshev-Fourier 

#General distribution generated by the product of Chebyshev and Fourier expansion terms
type ChebFourIdx <: RefractiveIndex
    #coeff_1 is a 1D array with the ith entry the coefficient for the T_{i-1}(r) term.
    coeff_1::Array{Float64,1}
    #Necessary condition for continuity at origin: coefficients for even orders of i in T_i(r)*cos(j*theta) or T_i(r)*sin(j*theta) must be 0, since the value of T_i(0) is nonzero, which would produce a discontinuity at the origin. So we define: coeff_c is a 2D array with the (i,j)th entry the coefficient for the T_{2*i-1}(r)*cos(j*theta) term, while coeff_s is a 2D array with the (i,j)th entry the coefficient for the T_{2*i-1}(r)*sin(j*theta) term
    coeff_c::Array{Float64,2}
    coeff_s::Array{Float64,2}
    
    #Construct by fieldname
    ChebFourIdx(;coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array(Float64,(0,0))) = new(coeffmatrix(coeff_1),coeffmatrix(coeff_c),coeffmatrix(coeff_s))
    #Construct by position
    ChebFourIdx(coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array(Float64,(0,0))) = new(coeffmatrix(coeff_1),coeffmatrix(coeff_c),coeffmatrix(coeff_s))
end

#Index value ------------------------------------------------------------------
function nfunc(idx::ChebFourIdx,r::Float64,theta::Float64)
    #Computing Chebyshev terms using recurrence relation
    n::Float64 = 0.0
    
    #Begin Tval at T_0
    Tval_2::Float64 = 2.0*sq(r)-1; Tval_1::Float64 = r; Tval::Float64=1.0
    
    for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
        #Compute Chebyshev value T_{i-1}(r) using recurrence relation
        Tval = 2.0*r*Tval_1 - Tval_2
        Tval_2 = Tval_1; Tval_1 = Tval
        
        if i <= length(idx.coeff_1)
            n += idx.coeff_1[i]*Tval
        end
        if mod(i,2) == 0
            for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                    n += idx.coeff_c[i/2,j]*Tval*cos(j*theta)
                end
                if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                    n += idx.coeff_s[i/2,j]*Tval*sin(j*theta)
                end
            end
        end
    end
    return n::Float64
end

function nfunc{N}(idx::ChebFourIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta))
    #Computing Chebyshev terms using recurrence relation
    n::Array{Float64,N} = zeros(Float64,size(r))
    
    #Begin Tval at T_0
    Tval_2::Array{Float64,N} = 2.0*sq(r)-1.0; Tval_1::Array{Float64,N} = r; Tval::Array{Float64,N}=ones(Float64,size(r))
    
    for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
        #Compute Chebyshev value T_{i-1}(r) using recurrence relation
        Tval = 2.0*r.*Tval_1 - Tval_2
        Tval_2 = Tval_1; Tval_1 = Tval
        
        if i <= length(idx.coeff_1)
            n += idx.coeff_1[i]*Tval
        end
        if mod(i,2) == 0
            for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                    n += idx.coeff_c[i/2,j]*Tval.*cos(j*theta)
                end
                if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                    n += idx.coeff_s[i/2,j]*Tval.*sin(j*theta)
                end
            end
        end
    end
    return n::Array{Float64,N}
end
        
#Index value and derivatives --------------------------------------------------
function nderiv(idx::ChebFourIdx,r::Float64,theta::Float64)
    n::Float64 = 0.0
    dr_n::Float64 = 0.0
    dtheta_n::Float64 = 0.0
    
    #Begin Tval and Uval at T_0 and U_{-1} respectively
    Tval_2::Float64 = 2.0*sq(r)-1; Tval_1::Float64 = r; Tval::Float64 = 1.0
    Uval_2::Float64 = -2.0*r; Uval_1::Float64 = -1.0; Uval::Float64 = 0.0
    
    for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
        #Compute Chebyshev value T_{i-1}(r) using recurrence relation
        Tval = 2.0*r*Tval_1 - Tval_2
        Tval_2 = Tval_1; Tval_1 = Tval
        #Compute (second kind) Chebyshev value U_{i-1}(r) using recurrence relation, to compute d(T_i)/dx = i*U_{i-1}
        Uval = 2.0*r*Uval_1 - Uval_2
        Uval_2 = Uval_1; Uval_1 = Uval
        
        if i <= length(idx.coeff_1)
            n += idx.coeff_1[i]*Tval
            dr_n += idx.coeff_1[i]*(i-1)*Uval; #Uval is one term behind
        end
        if mod(i,2) == 0
            for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                    n += idx.coeff_c[i/2,j]*Tval*cos(j*theta)
                    dr_n += idx.coeff_c[i/2,j]*(i-1)*Uval*cos(j*theta)
                    dtheta_n -= idx.coeff_c[i/2,j]*Tval*j*sin(j*theta)
                end
                if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                    n += idx.coeff_s[i/2,j]*Tval*sin(j*theta)
                    dr_n += idx.coeff_s[i/2,j]*(i-1)*Uval*sin(j*theta)
                    dtheta_n += idx.coeff_s[i/2,j]*Tval*j*cos(j*theta)
                end
            end
        end
    end
    return (n,dr_n,dtheta_n)
end

function nderiv{N}(idx::ChebFourIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta))
    n::Array{Float64,N} = similar(r)
    dr_n::Array{Float64,N} = similar(r)
    dtheta_n::Array{Float64,N} = similar(r)
    
    for k=1:length(r)
        n[k] = 0.0
        dr_n[k] = 0.0
        dtheta_n[k] = 0.0
        #Begin Tval and Uval at T_0 and U_{-1} respectively
        Tval_2::Float64 = 2.0*sq(r[k])-1.0; Tval_1::Float64 = r[k]; Tval::Float64 = 1.0
        Uval_2::Float64 = -2.0*r[k]; Uval_1::Float64 = -1.0; Uval::Float64 = 0.0
        
        for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
            #Compute Chebyshev value T_{i-1}(r) using recurrence relation
            Tval = 2.0*r[k]*Tval_1 - Tval_2
            Tval_2 = Tval_1; Tval_1 = Tval
            #Compute (second kind) Chebyshev value U_{i-1}(r) using recurrence relation, to compute d(T_i)/dx = i*U_{i-1}
            Uval = 2.0*r[k]*Uval_1 - Uval_2
            Uval_2 = Uval_1; Uval_1 = Uval
            
            if i <= length(idx.coeff_1)
                n[k] += idx.coeff_1[i]*Tval
                dr_n[k] += idx.coeff_1[i]*(i-1)*Uval; #Uval is one term behind
            end
            if mod(i,2) == 0
                for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                    if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                        n[k] += idx.coeff_c[i/2,j]*Tval*cos(j*theta[k])
                        dr_n[k] += idx.coeff_c[i/2,j]*(i-1)*Uval*cos(j*theta[k])
                        dtheta_n[k] -= idx.coeff_c[i/2,j]*j*Tval*sin(j*theta[k])
                    end
                    if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                        n[k] += idx.coeff_s[i/2,j]*Tval*sin(j*theta[k])
                        dr_n[k] += idx.coeff_s[i/2,j]*(i-1)*Uval*sin(j*theta[k])
                        dtheta_n[k] += idx.coeff_s[i/2,j]*j*Tval*cos(j*theta[k])
                    end
                end
            end
        end
    end
    return (n,dr_n,dtheta_n)
end


# #############################################################################
# #############################################################################
#Combinator refractive index distribution

#Sum of refractive index distributions
type SumIdx <: RefractiveIndex
    #List of refractive index distributions
    idx_list::Array{RefractiveIndex,1}
    
    SumIdx(idxs...) = new([idxs...])
end

function nfunc(idx::SumIdx,r::Float64,theta::Float64)
    n::Float64 = 0.0
    for i=1:length(idx.idx_list)
        n += nfunc(idx.idx_list[i],r,theta)
    end
    return n
end

function nderiv(idx::SumIdx,r::Float64,theta::Float64)
    N::Float64 = 0.0
    dr_N::Float64 = 0.0
    dtheta_N::Float64 = 0.0
    for i=1:length(idx.idx_list)
        n::Float64,dr_n::Float64,dtheta_n::Float64 = nderiv(idx.idx_list[i],r,theta)
        N += n
        dr_N += dr_n
        dtheta_N += dtheta_n
    end
    return (N,dr_N,dtheta_N)
end

