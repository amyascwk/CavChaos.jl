
#==============================================================================
#Helpful functions

sq{T<:Number}(x::T) = x*x::T;
cb{T<:Number}(x::T) = x*x*x::T;
@vectorize_1arg Number sq;
@vectorize_1arg Number cb;


#==============================================================================
#List of Boundaries


#Each boundary bnd has at theta a radial position of r = bnd.rfunc(theta) and an angle of normang between the normal vector and the horizontal. bnd.sys(theta) returns an array containing r and normang, while bnd.sys!(theta,A) stores the result in array A.

abstract Boundary;


#Circle
type CircularBnd <: Boundary
    name::String;
    rfunc::Function;
    sys::Function;
    sys!::Function;
    
    function CircularBnd()
        bnd = new("circular");
        
        #Radial position
        rfunc(theta::Float64) = 1.0::Float64;
        @vectorize_1arg Float64 rfunc;
        bnd.rfunc = rfunc;
        
        #Radial position and norm vector angle
        sys(theta::Float64) = [1.0,theta]::Array{Float64,1};
        function sys{N}(theta::Array{Float64,N})
            result::Array{Array{Float64,N}} = Array(Array{Float64,N},2);
            result[1] = ones(Float64,size(theta));
            result[2] = theta;
            return result;
        end
        bnd.sys = sys;
        
        function sys!(theta::Float64,result::Array{Float64,1})
            result[1],result[2] = bnd.sys(theta);
            return nothing;
        end
        bnd.sys! = sys!;
        
        return bnd;
    end
end


#Ellipse
type EllipticalBnd <: Boundary
    name::String;
    eccentricity::Float64;
    rfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    EllipticalBnd(;eccentricity::Float64=0.5) = EllipticalBnd(eccentricity);
    #Construct by abbreviation
    EllipticalBnd(;ecc::Float64=0.5) = EllipticalBnd(ecc);
    #Construct by position
    function EllipticalBnd(ecc::Float64=0.5)
        bnd = new("elliptical",ecc);
        
        #Radial position
        function rfunc(theta::Float64)
            #r = 1/sqrt(cos(theta)^2 + sin(theta)^2/(1-e^2))
            r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
            return r::Float64;
        end
        @vectorize_1arg Float64 rfunc;
        bnd.rfunc = rfunc;
        
        #Radial position and norm vector angle
        function sys(theta::Float64)
            r::Float64 = 1.0/sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
            dtheta_r::Float64 = -0.5*sin(2*theta)*sq(bnd.eccentricity)*cb(r)/(1-sq(bnd.eccentricity));
            normang::Float64 = theta - atan2(dtheta_r,r);
            return [r,normang]::Array{Float64,1};
        end
        function sys{N}(theta::Array{Float64,N})
            r::Array{Float64,N} = 1.0 ./sqrt(sq(cos(theta))+sq(sin(theta))/(1-sq(bnd.eccentricity)));
            dtheta_r::Array{Float64,N} = -0.5*sin(2*theta)*sq(bnd.eccentricity).*cb(r)/(1-sq(bnd.eccentricity));
            normang::Array{Float64,N} = theta - atan2(dtheta_r,r);
            
            result::Array{Array{Float64,N},1} = Array(Array{Float64,N},2);
            result[1] = r;
            result[2] = normang;
            return result;
        end
        bnd.sys = sys;
        
        function sys!(theta::Float64,result::Array{Float64,1})
            result[1],result[2] = bnd.sys(theta);
            return nothing;
        end
        bnd.sys! = sys!;
        
        return bnd;
    end
end


#Regular Polygon
type PolygonalBnd <: Boundary
    #roundedness is the degree of deformation into a circle
    name::String
    sides::Int64;
    roundedness::Float64;
    rfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    PolygonalBnd(;sides::Int64=3,roundedness::Float64=0.0) = PolygonalBnd(sides,roundedness);
    #Construct by abbreviation
    PolygonalBnd(;N::Int64=3,rnd::Float64=0.0) = PolygonalBnd(N,rnd);
    #Construct by position
    function PolygonalBnd(N::Int64=3,rnd::Float64=0.0)
        bnd = new("polygonal",N,rnd);
        
        #Radial position
        function rfunc(theta::Float64)
            angdiv::Float64 = 2*pi/bnd.sides;
            thetap::Float64 = mod(theta,angdiv)-0.5*angdiv;
            r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv);
            return r::Float64;
        end
        @vectorize_1arg Float64 rfunc;
        bnd.rfunc = rfunc;
        
        #Radial position and norm vector angle
        function sys(theta::Float64)
            angdiv::Float64 = 2*pi/bnd.sides;
            thetap::Float64 = mod(theta,angdiv)-0.5*angdiv;
            r::Float64 = (1-bnd.roundedness)/cos(thetap) + bnd.roundedness/cos(0.5*angdiv);
            normang::Float64 = theta;
            if mod(theta,angdiv) != 0.0
                dtheta_r::Float64 = (1-bnd.roundedness)*sin(thetap)/sq(cos(thetap));
                normang -= atan2(dtheta_r,r);
            end
            return [r,normang]::Array{Float64,1};
        end
        function sys{N}(theta::Array{Float64,N})
            angdiv::Array{Float64,N} = 2*pi/bnd.sides*ones(Float64,size(theta));
            thetap::Array{Float64,N} = mod(theta,angdiv)-0.5*angdiv;
            r::Array{Float64,N} = (1-bnd.roundedness)./cos(thetap) + bnd.roundedness./cos(0.5*angdiv);
            cornermiss::Array{Bool,N} = (mod(theta,angdiv) .!= 0.0);
            normang::Array{Float64,N} = theta;
            dtheta_r::Array{Float64,N} = (1-bnd.roundedness)*sin(thetap[cornermiss])./sq(cos(thetap[cornermiss]));
            normang[cornermiss] -= atan2(dtheta_r,r);
            
            result::Array{Array{Float64,N},1} = Array(Array{Float64,N},2);
            result[1] = r;
            result[2] = normang;
            return result;
        end
        bnd.sys = sys;
        
        function sys!(theta::Float64,result::Array{Float64,1})
            result[1],result[2] = bnd.sys(theta);
            return nothing;
        end
        bnd.sys! = sys!;
        
        return bnd;
    end
end


#General boundary from Fourier expansion
type FourierBnd <: Boundary
    #coeff_c and coeff_s are vectors respectively containing the Fourier coefficients for the cosine and sine series expansion of the radius function, normalized to have a constant component of 1.0
    name::String;
    coeff_c::Array{Float64,1};
    coeff_s::Array{Float64,1};
    rfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    FourierBnd(;coeff_c::Array{Float64,1}=[0.0,0.02],coeff_s::Array{Float64,1}=Float64[]) = GeneralBnd(coeff_c,coeff_s);
    #Construct by position
    function FourierBnd(coeff_c::Array{Float64,1}=[0.0,0.05],coeff_s::Array{Float64,1}=Float64[])
        bnd = new("fourier",coeff_c,coeff_s);
        
        #Radial position
        function rfunc(theta::Float64)
            r::Float64 = 1.0;
            for i = 1:length(bnd.coeff_c)
                r += bnd.coeff_c[i]*cos(i*theta);
            end
            for i = 1:length(bnd.coeff_s)
                r += bnd.coeff_s[i]*sin(i*theta);
            end
            return r::Float64;
        end
        @vectorize_1arg Float64 rfunc;
        bnd.rfunc = rfunc;
        
        #Radial position and norm vector angle
        function sys(theta::Float64)
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
            return [r,normang]::Array{Float64,1};
        end
        function sys{N}(theta::Array{Float64,N})
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
            result::Array{Array{Float64,N},1} = Array(Array{Float64,N},2);
            result[1] = r;
            result[2] = normang;
            return result::Array{Array{Float64,N},1};
        end
        bnd.sys = sys;
        
        function sys!(theta::Float64,result::Array{Float64,1})
            result[1],result[2] = bnd.sys(theta);
            return nothing;
        end
        bnd.sys! = sys!;
        
        return bnd;
    end
end
