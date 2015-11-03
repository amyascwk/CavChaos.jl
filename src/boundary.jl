# #############################################################################
# #############################################################################
#Cavity boundary shape

#This file contains functions for specifying the shape of the cavity boundary, by subtypes of the abstract type Boundary. The following are the available constructors, which can take both positional or keyword arguments:

#   CircularBnd()
#   constructs a circular boundary of radius 1.0

#   EllipticalBnd(eccentricity::Float64=0.5)
#   constructs an elliptical boundary of semimajor axis 1.0 and a specified eccentricity (default 0.5).

#   PolygonalBnd(sides::Int64=3,roundness::Float64=0.0)
#   constructs a regular polygonal boundary of a specified number of sides (default 3) and amount of deformation into a circle (0.0 to 1.0, default 0.0).

#   FourierBnd(coeff_c::Array{Float64,1}=[0.0,0.05], coeff_s::Array{Float64,1}=Float64[])
#   constructs a general boundary based on the specified Fourier expansion coefficients of the boundary radial position in polar coordinates, with a constant component of 1.0 (default is r = 1.0 + 0.05*cos(2*theta)).

#Information about the boundary in polar coordinates can be extracted by the following functions, which dispatch according to the boundary subtype to the appropriate calculations for the shape:

#   r = rfunc(bnd::Boundary,theta)
#   returns the radial position of the boundary <bnd> at the angular position <theta>. <theta> can be a Float64 or Array{Float64,N}, and the result will be of the same type and dimensions as <theta>.

#   r,dtheta_r = rderiv(bnd::Boundary,theta)
#   returns a tuple (<r>,<dtheta_r>) of the radial position <r> and its derivative <dtheta_r> with respect to angular position <theta>, for the boundary <bnd> at the angular position <theta>. <theta> can be a Float64 or Array{Float64,N}, and each entry in the tuple result will be of the same type and dimensions as <theta>. There are no functions to compute the derivatives alone since the calculations involved almost always requires the radial position to also be computed.

#   r,dtheta_r,dtheta2_r = rderiv2(bnd::Boundary,theta)
#   this is like rderiv but also includes the second derivative <dtheta2_r> in the resulting tuple.

#   alpha = normang(theta,r,dtheta_r)
#   returns the angle <alpha> of the normal vector at angular position <theta>. <theta> can be a Float64 or Array{Float64,N}, but <r> and <dtheta_r> must follow its type and dimensions, and the result <alpha> will also follow its type and dimensions.

#   r,alpha = rsys(bnd::Boundary,theta)
#   returns the tuple (<r>,<alpha>) of the radial position r and the angle alpha of the normal vector, for the boundary bnd at the angular position theta. Theta can be a Float64 or Array{Float64,N}, and each entry in the tuple result will be of the same type and dimensions as theta. This function is provided because the radial position and normal vector angle are the particular quantities of interest.

#   rsys!(bnd:Boundary,theta::Float64,result::Array{Float64,1})
#   This is the in-place version of rsys storing the results as the first two entries of the array <result>.


# #############################################################################
# #############################################################################
#Initialize

#Dependencies
#require("util.jl")
import Base: ==, hash

#General boundary type
abstract Boundary


# #############################################################################
# #############################################################################
#Hashing and equality

#Extend base hash function with method to hash equivalent Boundary objects identically.
hash(bnd::Boundary) = fieldhashing(bnd)
hash(bnd::Boundary,h::UInt64) = hash(hash(bnd),h)

#Equality is equality of type and field values
==(bnd1::Boundary,bnd2::Boundary) = fieldequality(bnd1,bnd2)


# #############################################################################
# #############################################################################
#Functions for generic boundary

#These functions are vectorized versions of rfunc, rderiv and rderiv2 for a general Boundary, as well as (vectorized/in-place versions of) rsys for use in solving for the ray trajectories. The implementations for specific Boundary subtypes, if present, are slightly more efficient.

#Radial position --------------------------------------------------------------
function rfunc{N}(bnd::Boundary,theta::Array{Float64,N})
    #Initialize array
    r::Array{Float64,N} = similar(theta)
    #Set values to rfunc output
    for i=1:length(theta)
        r[i] = rfunc(bnd,theta[i])
    end
    return r::Array{Float64,N}
end

#Radial position and derivatives ----------------------------------------------
function rderiv{N}(bnd::Boundary,theta::Array{Float64,N})
    #Initialize arrays
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    #Set values to rderiv output
    for i=1:length(theta)
        r[i],dtheta_r[i] = rderiv(bnd,theta[i])
    end
    return (r,dtheta_r)
end

#Radial position and derivatives2 ---------------------------------------------
function rderiv2{N}(bnd::Boundary,theta::Array{Float64,N})
    #Initialize arrays
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    dtheta2_r::Array{Float64,N} = similar(theta)
    #Set values to rderiv2 output
    for i=1:length(theta)
        r[i],dtheta_r[i],dtheta2_r[i] = rderiv2(bnd,theta[i])
    end
    return (r,dtheta_r,dtheta2_r)
end


#Angular position of normal vector --------------------------------------------
function normang(theta::Float64,r::Float64,dtheta_r::Float64)
    return (theta - atan2(dtheta_r,r))::Float64
end

function normang{N}(theta::Array{Float64,N},r::Array{Float64,N},dtheta_r::Array{Float64,N})
    #Initialize arrays
    alpha::Array{Float64,N}  = similar(theta)
    #Set values to normang output
    for i=1:length(theta)
        alpha[i] = normang(theta[i],r[i],dtheta_r[i])
    end
    return alpha
end


#Radial position and normal vector angle --------------------------------------
function rsys(bnd::Boundary,theta::Float64)
    r::Float64, dtheta_r::Float64 = rderiv(bnd,theta)
    alpha::Float64 = normang(theta,r,dtheta_r)
    return (r,alpha)
end

function rsys{N}(bnd::Boundary,theta::Array{Float64,N})
    #Initialize arrays
    r::Array{Float64,N} = similar(theta)
    alpha::Array{Float64,N} = similar(theta)
    #Set values to rsys output
    for i=1:length(theta)
        r[i],alpha[i] = rsys(bnd,theta[i])
    end
    return (r,alpha)
end

#In-place storing version of rsys
function rsys!(bnd::Boundary,theta::Float64,result::Array{Float64,1})
    result[1],dtheta_r::Float64 = rderiv(bnd,theta)
    result[2] = normang(theta,result[1],dtheta_r)
    return nothing
end


# #############################################################################
# #############################################################################
#List of Boundaries

#Code is replicated a lot here to reduce the overhead of function calls, because these calculations are repeated a lot, so speed is fairly important.

# #############################################################################
#Circle

#Circular Boundary of radius 1.0
type CircularBnd <: Boundary
    #Contains no fields
end

#Radial position --------------------------------------------------------------
rfunc(bnd::CircularBnd,theta::Float64) = 1.0::Float64
rfunc{N}(bnd::CircularBnd,theta::Array{Float64,N}) = ones(Float64,size(theta))::Array{Float64,N}

#Radial position and derivatives ----------------------------------------------
rderiv(bnd::CircularBnd,theta::Float64) = (1.0,0.0)
function rderiv{N}(bnd::CircularBnd,theta::Array{Float64,N})
    return (ones(Float64,size(theta)),zeros(Float64,size(theta)))
end

#Radial position and derivatives2 ---------------------------------------------
rderiv2(bnd::CircularBnd,theta::Float64) = (1.0,0.0,0.0)
function rderiv2{N}(bnd::CircularBnd,theta::Array{Float64,N})
    return (ones(Float64,size(theta)),zeros(Float64,size(theta)),zeros(Float64,size(theta)))
end


# #############################################################################
#Ellipse

#Elliptical Boundary with semimajor axis of length 1.0
type EllipticalBnd <: Boundary
    #Eccentricity of the ellipse
    eccentricity::Float64
    
    #Construct by fieldname
    EllipticalBnd(;eccentricity::Float64=0.5) = new(eccentricity)
    #Construct by position
    EllipticalBnd(eccentricity::Float64=0.5) = new(eccentricity)
end

#Radial position --------------------------------------------------------------
function rfunc(bnd::EllipticalBnd,theta::Float64)
    r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)))
    return r::Float64
end

function rfunc{N}(bnd::EllipticalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    for i=1:length(theta)
        r[i] = 1.0/sqrt(sq(cos(theta[i]))+sq(sin(theta[i]))/(1-sq(bnd.eccentricity)))
    end
    return r::Array{Float64,N}
end

#Radial position and derivatives ----------------------------------------------
function rderiv(bnd::EllipticalBnd,theta::Float64)
    r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)))
    dtheta_r::Float64 = -0.5*sin(2*theta)*sq(bnd.eccentricity)*cb(r)/(1-sq(bnd.eccentricity))
    return (r,dtheta_r)
end

function rderiv{N}(bnd::EllipticalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    for i=1:length(theta)
        r[i] = 1.0/sqrt(sq(cos(theta[i]))+sq(sin(theta[i]))/(1-sq(bnd.eccentricity)))
        dtheta_r[i] = -0.5*sin(2*theta[i])*sq(bnd.eccentricity)*cb(r[i])/(1-sq(bnd.eccentricity))
    end
    return (r,dtheta_r)
end

#Radial position and derivatives2 ---------------------------------------------
function rderiv2(bnd::EllipticalBnd,theta::Float64)
    r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)))
    dtheta_r::Float64 = -0.5*sin(2*theta)*sq(bnd.eccentricity)*cb(r)/(1-sq(bnd.eccentricity))
    dtheta2_r::Float64 = -(cos(2*theta)*r + 1.5*sin(2*theta)*dtheta_r)*sq(r)*sq(bnd.eccentricity)/(1-sq(bnd.eccentricity))
    return (r,dtheta_r,dtheta2_r)
end

function rderiv2{N}(bnd::EllipticalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    dtheta2_r::Array{Float64,N} = similar(theta)
    for i=1:length(theta)
        r[i] = 1.0/sqrt(sq(cos(theta[i]))+sq(sin(theta[i]))/(1-sq(bnd.eccentricity)))
        dtheta_r[i] = -0.5*sin(2*theta[i])*sq(bnd.eccentricity)*cb(r[i])/(1-sq(bnd.eccentricity))
        dtheta2_r[i] = -(cos(2*theta[i])*r[i] + 1.5*sin(2*theta[i])*dtheta_r[i])*sq(r[i])*sq(bnd.eccentricity)/(1-sq(bnd.eccentricity))
    end
    return (r,dtheta_r,dtheta2_r)
end


# #############################################################################
#Polygon

#Regular Polygon with center to vertex distance of 1.0
type PolygonalBnd <: Boundary
    #Number of sides of the regular polygon
    sides::Int64
    #Roundness (the degree of deformation into a circle) taking a value of 0.0 for a polygon to 1.0 for a circle.
    roundedness::Float64; 
    
    #Construct by fieldname
    PolygonalBnd(;sides::Int64=3,roundness::Float64=0.0) = new(sides,roundness)
    #Construct by position
    PolygonalBnd(sides::Int64=3,roundness::Float64=0.0) = new(sides,roundness)
end
        
#Radial position --------------------------------------------------------------
function rfunc(bnd::PolygonalBnd,theta::Float64)
    angdiv::Float64 = 2*pi/bnd.sides
    thetap::Float64 = mod(theta,angdiv)-0.5*angdiv
    r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv)
    return r::Float64
end

function rfunc{N}(bnd::PolygonalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    angdiv::Float64 = 2*pi/bnd.sides
    for i=1:length(theta)
        thetap::Float64 = mod(theta[i],angdiv)-0.5*angdiv
        r[i] = (1-bnd.roundedness)./cos(thetap) + bnd.roundedness./cos(0.5*angdiv)
    end
    return r::Array{Float64,N}
end

#Radial position and derivatives ----------------------------------------------
function rderiv(bnd::PolygonalBnd,theta::Float64)
    angdiv::Float64 = 2*pi/bnd.sides
    thetap::Float64 = mod(theta,angdiv)-0.5*angdiv
    r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv)
    if mod(theta,angdiv) != 0.0
        dtheta_r::Float64 = (1-bnd.roundedness)*sec(thetap)*tan(thetap)
    else
        dtheta_r = 0.0
    end
    return (r,dtheta_r)
end

function rderiv{N}(bnd::PolygonalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    angdiv::Float64 = 2*pi/bnd.sides
    for i=1:length(theta)
        thetap::Float64 = mod(theta[i],angdiv)-0.5*angdiv
        r[i] = (1-bnd.roundedness)./cos(thetap) + bnd.roundedness./cos(0.5*angdiv)
        if mod(theta[i],angdiv) != 0.0
            dtheta_r[i] = (1-bnd.roundedness)*sec(thetap)*tan(thetap)
        else
            dtheta_r[i] = 0.0
        end
    end
    return (r,dtheta_r)
end

#Radial position and derivatives2 ---------------------------------------------
function rderiv2(bnd::PolygonalBnd,theta::Float64)
    angdiv::Float64 = 2*pi/bnd.sides
    thetap::Float64 = mod(theta,angdiv)-0.5*angdiv
    r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv)
    if mod(theta,angdiv) != 0.0
        dtheta_r::Float64 = (1-bnd.roundedness)*sec(thetap)*tan(thetap)
        dtheta2_r::Float64 = (1-bnd.roundedness)*sec(thetap)*(1+2*sq(tan(thetap)))
    else
        dtheta_r = 0.0
        dtheta2_r = 0.0
    end
    return (r,dtheta_r,dtheta2_r)
end

function rderiv2{N}(bnd::PolygonalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    dtheta2_r::Array{Float64,N} = similar(theta)
    angdiv::Float64 = 2*pi/bnd.sides
    for i=1:length(theta)
        thetap::Float64 = mod(theta[i],angdiv)-0.5*angdiv
        r[i] = (1-bnd.roundedness)./cos(thetap) + bnd.roundedness./cos(0.5*angdiv)
        if mod(theta[i],angdiv) != 0.0
            dtheta_r[i] = (1-bnd.roundedness)*sec(thetap)*tan(thetap)
            dtheta2_r[i] = (1-bnd.roundedness)*sec(thetap)*(1+2*sq(tan(thetap)))
        else
            dtheta_r[i] = 0.0
            dtheta2_r[i] = 0.0
        end
    end
    return (r,dtheta_r,dtheta2_r)
end


# #############################################################################
#Fourier

#General boundary generated by Fourier expansion, normalized to have a constant component of 1.0
type FourierBnd <: Boundary
    #coeff_c and coeff_s are arrays containing respectively the Fourier coefficients for the cosine and sine series expansion of the radius function. Specifically, the nth value of coeff_c is the coefficient of the cos(n*theta) term, etc. They need only be long enough to include the last nonzero coefficient of each series.
    coeff_c::Array{Float64,1}
    coeff_s::Array{Float64,1}
    
    #Construct by fieldname
    FourierBnd(;coeff_c::Array{Float64,1}=[0.0,0.05],coeff_s::Array{Float64,1}=Float64[]) = new(coeffmatrix(coeff_c),coeffmatrix(coeff_s))
    #Construct by position
    FourierBnd(coeff_c::Array{Float64,1}=[0.0,0.05],coeff_s::Array{Float64,1}=Float64[]) = new(coeffmatrix(coeff_c),coeffmatrix(coeff_s))
end

#Radial position --------------------------------------------------------------
function rfunc(bnd::FourierBnd,theta::Float64)
    r::Float64 = 1.0
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta)
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta)
    end
    return r::Float64
end

function rfunc{N}(bnd::FourierBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    for i=1:length(theta)
        r[i] = 1.0
        for j = 1:length(bnd.coeff_c)
            r[i] += bnd.coeff_c[j]*cos(j*theta[i])
        end
        for j = 1:length(bnd.coeff_s)
            r[i] += bnd.coeff_s[j]*sin(j*theta[i])
        end
    end
    return r::Array{Float64,N}
end

#Radial position and derivatives ----------------------------------------------
function rderiv(bnd::FourierBnd,theta::Float64)
    r::Float64 = 1.0
    dtheta_r::Float64 = 0.0
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta)
        dtheta_r -= bnd.coeff_c[i]*i*sin(i*theta)
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta)
        dtheta_r += bnd.coeff_s[i]*i*cos(i*theta)
    end
    return (r,dtheta_r)
end

function rderiv{N}(bnd::FourierBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    for i=1:length(theta)
        r[i] = 1.0
        dtheta_r[i] = 0.0
        for j = 1:length(bnd.coeff_c)
            r[i] += bnd.coeff_c[j]*cos(j*theta[i])
            dtheta_r[i] -= bnd.coeff_c[j]*j*sin(j*theta[i])
        end
        for j = 1:length(bnd.coeff_s)
            r[i] += bnd.coeff_s[j]*sin(j*theta[i])
            dtheta_r[i] += bnd.coeff_s[j]*j*cos(j*theta[i])
        end
    end
    return (r,dtheta_r)
end

#Radial position and derivatives2 ---------------------------------------------
function rderiv2(bnd::FourierBnd,theta::Float64)
    r::Float64 = 1.0
    dtheta_r::Float64 = 0.0
    dtheta2_r::Float64 = 0.0
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta)
        dtheta_r -= bnd.coeff_c[i]*i*sin(i*theta)
        dtheta2_r -= bnd.coeff_c[i]*sq(i)*cos(i*theta)
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta)
        dtheta_r += bnd.coeff_s[i]*i*cos(i*theta)
        dtheta2_r -= bnd.coeff_s[i]*sq(i)*sin(i*theta)
    end
    return (r,dtheta_r,dtheta2_r)
end

function rderiv2{N}(bnd::FourierBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = similar(theta)
    dtheta_r::Array{Float64,N} = similar(theta)
    dtheta2_r::Array{Float64,N} = similar(theta)
    for i=1:length(theta)
        r[i] = 1.0
        dtheta_r[i] = 0.0
        dtheta2_r[i] = 0.0
        for j = 1:length(bnd.coeff_c)
            r[i] += bnd.coeff_c[j]*cos(j*theta[i])
            dtheta_r[i] -= bnd.coeff_c[j]*j*sin(j*theta[i])
            dtheta2_r[i] -= bnd.coeff_c[j]*sq(j)*cos(j*theta[i])
        end
        for j = 1:length(bnd.coeff_s)
            r[i] += bnd.coeff_s[j]*sin(j*theta[i])
            dtheta_r[i] += bnd.coeff_s[j]*j*cos(j*theta[i])
            dtheta2_r[i] -= bnd.coeff_s[j]*sq(j)*sin(j*theta[i])
        end
    end
    return (r,dtheta_r,dtheta2_r)
end
