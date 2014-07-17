

#==============================================================================
#Helpful functions

sq{T<:Number}(x::T) = x*x::T;
cb{T<:Number}(x::T) = x*x*x::T;
@vectorize_1arg Number sq;
@vectorize_1arg Number cb;


#==============================================================================
#List of Boundaries


#Each boundary bnd has at theta a radial position of r = rfunc(bnd,theta) and an angle of normang between the normal vector and the horizontal. rsys(bnd,theta) returns an array containing r and normang, while rsys!(bnd,theta,A) stores the result in array A.

abstract Boundary;

function rsys!(bnd::Boundary,theta::Float64,result::Array{Float64,1})
    result[1],result[2] = rsys(bnd,theta);
    return nothing;
end


#Circle -----------------------------------------------------------------------
type CircularBnd <: Boundary
    #No arguments
end

#Radial position
rfunc(bnd::CircularBnd,theta::Float64) = 1.0::Float64;
rfunc{N}(bnd::CircularBnd,theta::Array{Float64,N}) = ones(Float64,size(theta))::Array{Float64,N};

#Radial position and norm vector angle
rsys(bnd::CircularBnd,theta::Float64) = (1.0,theta);
function rsys{N}(bnd::CircularBnd,theta::Array{Float64,N})
    return (ones(Float64,size(theta)),copy(theta));
end


#Ellipse ----------------------------------------------------------------------
type EllipticalBnd <: Boundary
    #1 argument: eccentricity of the ellipse
    eccentricity::Float64;
    
    #Construct by fieldname
    EllipticalBnd(;eccentricity::Float64=0.5) = new(eccentricity);
    #Construct by abbreviation
    EllipticalBnd(;ecc::Float64=0.5) = new(ecc);
    #Construct by position
    EllipticalBnd(ecc::Float64=0.5) = new(ecc);
end

#Radial position
#r = 1/sqrt(cos(theta)^2 + sin(theta)^2/(1-e^2))
function rfunc(bnd::EllipticalBnd,theta::Float64)
    r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
    return r::Float64;
end

function rfunc{N}(bnd::EllipticalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = 1.0./sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
    return r::Array{Float64,N};
end
        
#Radial position and norm vector angle
function rsys(bnd::EllipticalBnd,theta::Float64)
    r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
    dtheta_r::Float64 = -0.5*sin(2*theta)*sq(bnd.eccentricity)*cb(r)/(1-sq(bnd.eccentricity));
    normang::Float64 = theta - atan2(dtheta_r,r);
    return (r,normang);
end

function rsys{N}(bnd::EllipticalBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = 1.0 ./sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
    dtheta_r::Array{Float64,N} = -0.5*sin(2*theta)*sq(bnd.eccentricity).*cb(r)/(1-sq(bnd.eccentricity));
    normang::Array{Float64,N} = theta - atan2(dtheta_r,r);
    return (r,normang);
end


#Regular Polygon --------------------------------------------------------------
type PolygonalBnd <: Boundary
    #2 arguments: number of sides of the regular polygon, and roundness (the degree of 
    #deformation into a circle) taking a value of 0.0 (polygon) to 1.0 (circle).
    sides::Int64;
    roundedness::Float64; 
    
    #Construct by fieldname
    PolygonalBnd(;sides::Int64=3,roundedness::Float64=0.0) = new(sides,roundness);
    #Construct by abbreviation
    PolygonalBnd(;N::Int64=3,rnd::Float64=0.0) = new(N,rnd);
    #Construct by position
    PolygonalBnd(N::Int64=3,rnd::Float64=0.0) = new(N,rnd);
end
        
#Radial position
function rfunc(bnd::PolygonalBnd,theta::Float64)
    angdiv::Float64 = 2*pi/bnd.sides;
    thetap::Float64 = mod(theta,angdiv)-0.5*angdiv;
    r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv);
    return r::Float64;
end

function rfunc{N}(bnd::PolygonalBnd,theta::Array{Float64,N})
    angdiv::Array{Float64,N} = fill(2*pi/bnd.sides,size(theta));
    thetap::Array{Float64,N} = mod(theta,angdiv)-0.5*angdiv;
    r::Array{Float64,N} = (1-bnd.roundedness)./cos(thetap) + bnd.roundedness./cos(0.5*angdiv);
    return r::Array{Float64,N};
end

#Radial position and norm vector angle
function rsys(bnd::PolygonalBnd,theta::Float64)
    angdiv::Float64 = 2*pi/bnd.sides;
    thetap::Float64 = mod(theta,angdiv)-0.5*angdiv;
    r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv);
    normang::Float64 = theta;
    if mod(theta,angdiv) != 0.0
        dtheta_r::Float64 = (1-bnd.roundedness)*sin(thetap)/sq(cos(thetap));
        normang -= atan2(dtheta_r,r);
    end
    return (r,normang);
end

function rsys{N}(bnd::PolygonalBnd,theta::Array{Float64,N})
    angdiv::Array{Float64,N} = fill(2*pi/bnd.sides,size(theta));
    thetap::Array{Float64,N} = mod(theta,angdiv)-0.5*angdiv;
    r::Array{Float64,N} = (1-bnd.roundedness)./cos(thetap) + bnd.roundedness./cos(0.5*angdiv);
    cornermiss::Array{Bool,N} = (mod(theta,angdiv) .!= 0.0);
    normang::Array{Float64,N} = copy(theta);
    dtheta_r::Array{Float64,1} = (1-bnd.roundedness)*sin(thetap[cornermiss])./sq(cos(thetap[cornermiss]));
    normang[cornermiss] -= atan2(dtheta_r,r[cornermiss]);
    return (r,normang);
end


#General boundary from Fourier expansion --------------------------------------
type FourierBnd <: Boundary
    #2 arguments: coeff_c and coeff_s are vectors respectively containing the Fourier coefficients for the cosine and sine series expansion of the radius function, normalized to have a constant component of 1.0
    coeff_c::Array{Float64,1};
    coeff_s::Array{Float64,1};
    
    #Construct by fieldname
    FourierBnd(;coeff_c::Array{Float64,1}=[0.0,0.05],coeff_s::Array{Float64,1}=Float64[]) = new(coeff_c,coeff_s);
    #Construct by abbreviation
    FourierBnd(;cc::Array{Float64,1}=[0.0,0.05],cs::Array{Float64,1}=Float64[]) = new(cc,cs);
    #Construct by position
    FourierBnd(coeff_c::Array{Float64,1}=[0.0,0.05],coeff_s::Array{Float64,1}=Float64[]) = new(coeff_c,coeff_s);
end

#Radial position
function rfunc(bnd::FourierBnd,theta::Float64)
    r::Float64 = 1.0;
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta);
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta);
    end
    return r::Float64;
end

function rfunc{N}(bnd::FourierBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = ones(Float64,size(theta));
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta);
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta);
    end
    return r::Array{Float64,N};
end

#Radial position and norm vector angle
function rsys(bnd::FourierBnd,theta::Float64)
    r::Float64 = 1.0;
    dtheta_r::Float64 = 0.0;
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta);
        dtheta_r -= bnd.coeff_c[i]*i*sin(i*theta);
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta);
        dtheta_r += bnd.coeff_s[i]*i*cos(i*theta);
    end
    normang::Float64 = theta - atan2(dtheta_r,r);
    return (r,normang);
end

function rsys{N}(bnd::FourierBnd,theta::Array{Float64,N})
    r::Array{Float64,N} = ones(Float64,size(theta));
    dtheta_r::Array{Float64,N} = zeros(Float64,size(theta));
    for i = 1:length(bnd.coeff_c)
        r += bnd.coeff_c[i]*cos(i*theta);
        dtheta_r -= bnd.coeff_c[i]*i*sin(i*theta);
    end
    for i = 1:length(bnd.coeff_s)
        r += bnd.coeff_s[i]*sin(i*theta);
        dtheta_r += bnd.coeff_s[i]*i*cos(i*theta);
    end
    normang::Array{Float64,N} = theta - atan2(dtheta_r,r);
    return (r,normang);
end


