require("cavbnds.jl");
import GSL;

#====================================================================
#List of refractive index distributions


#Each refractive index distribution idx has at position (r,theta) a refractive index n = idx.nfunc(r,theta), with derivatives dr_n and dtheta_n for the radial and angular components respectively. idx.sys(r,theta) returns an array containing n, dr_n and dtheta_n, while idx.sys!(r,theta,A) stores the result in array A.

abstract RefractiveIndex;


#Uniform
type UniformIdx <: RefractiveIndex
    name::String;
    index_base::Float64
    nfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    UniformIdx(;index_base::Float64=2.0) = UniformIdx(index_base);
    #Construct by abbreviation
    UniformIdx(n0::Float64=2.0) = UniformIdx(n0);
    #Construct by position
    function UniformIdx(n0::Float64=2.0)
        idx = new("uniform",n0);
        
        #Index value
        nfunc(r::Float64,theta::Float64) = idx.index_base::Float64;
        @vectorize_2arg Float64 nfunc;
        idx.nfunc = nfunc;
        
        #Index value and derivatives
        sys(r::Float64,theta::Float64) = [idx.index_base,0.0,0.0]::Array{Float64,1};
        @vectorize_2arg Float64 sys;
        idx.sys = sys;
        
        function sys!(r::Float64,theta::Float64,result::Array{Float64,1})
            result[1],result[2],result[3] = idx.sys(r,theta);
            return nothing;
        end
        idx.sys! = sys!;
        
        return idx;
    end
end


#Gaussian
type GaussianIdx <: RefractiveIndex
    #The index_variation, spread and center fields are arrays where each row corresponds to the center positions, amplitude and amount of spread of a Gaussian variation contribution. Each row of the spread has 3 components: 2 for the standard deviations in 2 perpendicular directions and the last for the rotation angle.
    name::String;
    index_base::Float64;
    center::Array{Float64,2};
    index_variation::Array{Float64,1};
    spread::Array{Float64,2};
    nfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    GaussianIdx(;index_base::Float64=2.0,center::Array{Float64,2}=[0.0 0.0],index_variation::Array{Float64,1}=[-0.2],spread::Array{Float64,2}=[0.5 0.5 0.0]) = GaussianIdx(index_base,center,index_variation,spread);
    #Construct by abbreviation
    GaussianIdx(;n0::Float64=2.0,cpos::Array{Float64,2}=[0.0 0.0],dn::Array{Float64,1}=[-0.2],sig::Array{Float64,2}=[0.5 0.5 0.0]) = GaussianIdx(n0,cpos,dn,sig);
    #Construct by position
    function GaussianIdx(n0::Float64=2.0,cpos::Array{Float64,2}=[0.0 0.0],dn::Array{Float64,1}=[-0.2],sig::Array{Float64,2}=[0.5 0.5 0.0])
        
        if !(length(dn) == size(sig,1) == size(cpos,1))
            error("Index variation, spread and center arrays must be of the same number of rows.");
        end
        
        idx = new("gaussian",n0,cpos,dn,sig);
        
        #Index value
        function nfunc(r::Float64,theta::Float64)
            #assert (length(idx.index_variation) == size(idx.center,1) == size(idx.spread,1));
            x::Float64 = r*cos(theta);
            y::Float64 = r*sin(theta);
            n::Float64 = idx.index_base;
            
            for i = 1:length(idx.index_variation)
                dxr::Float64 = (x-idx.center[i,1])*cos(idx.spread[i,3]) - (y-idx.center[i,2])*sin(idx.spread[i,3]);
                dyr::Float64 = (y-idx.center[i,2])*cos(idx.spread[i,3]) + (x-idx.center[i,1])*sin(idx.spread[i,3]);
                n += idx.index_variation[i]*exp(-0.5*(sq(dxr/idx.spread[i,1]) + sq(dyr/idx.spread[i,2])));
            end
            return n;
        end
        @vectorize_2arg Float64 nfunc;
        idx.nfunc = nfunc;
        
        #Index value and derivatives
        function sys(r::Float64,theta::Float64)
            #assert (length(idx.index_variation) == size(idx.center,1) == size(idx.spread,1));
            x::Float64 = r*cos(theta);
            y::Float64 = r*sin(theta);
            n::Float64 = idx.index_base;
            dr_n::Float64 = 0.0;
            dtheta_n::Float64 = 0.0;
            
            for i = 1:length(idx.index_variation)
                dxr::Float64 = (x-idx.center[i,1])*cos(idx.spread[i,3]) - (y-idx.center[i,2])*sin(idx.spread[i,3]);
                dyr::Float64 = (y-idx.center[i,2])*cos(idx.spread[i,3]) + (x-idx.center[i,1])*sin(idx.spread[i,3]);
                gaussians::Float64 = idx.index_variation[i]*exp(-0.5*(sq(dxr/idx.spread[i,1]) + sq(dyr/idx.spread[i,2])));
                
                n += gaussians;
                dr_n -= (cos(theta+idx.spread[i,3])*dxr/sq(idx.spread[i,1]) + sin(theta+idx.spread[i,3])*dyr/sq(idx.spread[i,2]))*gaussians;
                dtheta_n += ((y*cos(idx.spread[i,3])+x*sin(idx.spread[i,3]))*dxr/sq(idx.spread[i,1]) + (-x*cos(idx.spread[i,3])+y*sin(idx.spread[i,3]))*dyr/sq(idx.spread[i,2]))*gaussians;
            end
            return [n,dr_n,dtheta_n]::Array{Float64,1};
        end
        #@vectorize_2arg Float64 sys;
        idx.sys = sys;
        
        function sys!(r::Float64,theta::Float64,result::Array{Float64,1})
            result[1],result[2],result[3] = idx.sys(r,theta);
            return nothing;
        end
        idx.sys! = sys!;
        
        return idx;
    end
end


#General index variation from Taylor expansion
type TaylorIdx <: RefractiveIndex
    # ***************** TO DO: Use sparse coefficients matrix? ****************
    #coeff is a matrix where the (i,j)th element is the coefficient of the x^i*y^j 
    #term in the Taylor expansion about the origin of the index distribution
    name::String;
    coeff::Array{Float64,2};
    nfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    TaylorIdx(;coeff::Array{Float64,2}=[2.0]) = TaylorIdx(coeff);
    #Construct by position
    function TaylorIdx(coeff::Array{Float64,2}=
        [2.0    0.0    0.1;
         0.0    0.0    0.01;
         0.1    0.01   0.1])
        
        idx = new("taylor",coeff);
        
        #Index value
        function nfunc(r::Float64,theta::Float64)
            x::Float64 = r*cos(theta);
            y::Float64 = r*sin(theta);
            n::Float64 = 0.0;
            for i = 1:size(idx.coeff,1)
                for j = 1:size(idx.coeff,2)
                    n += idx.coeff[i,j]*GSL.pow_int(y,i-1)*GSL.pow_int(x,j-1);
                end
            end
            return n;
        end
        @vectorize_2arg Float64 nfunc;
        idx.nfunc = nfunc;
        
        #Index value and derivatives
        function sys(r::Float64,theta::Float64)
            x::Float64 = r*cos(theta);
            y::Float64 = r*sin(theta);
            n::Float64 = 0.0;
            dr_n::Float64 = 0.0;
            dtheta_n::Float64 = 0.0;
            
            for i = 1:size(idx.coeff,1)
                for j = 1:size(idx.coeff,2)
                    n += idx.coeff[i,j]*GSL.pow_int(y,i-1)*GSL.pow_int(x,j-1);
                    if i >= 2
                        dr_n += idx.coeff[i,j]*(i-1)*GSL.pow_int(y,i-2)*GSL.pow_int(x,j-1)*sin(theta);
                        dtheta_n += idx.coeff[i,j]*(i-1)*GSL.pow_int(y,i-2)*GSL.pow_int(x,j);
                    end
                    if j >= 2
                        dr_n += idx.coeff[i,j]*GSL.pow_int(y,i-1)*(j-1)*GSL.pow_int(x,j-2)*cos(theta);
                        dtheta_n -= idx.coeff[i,j]*GSL.pow_int(y,i)*(j-1)*GSL.pow_int(x,j-2);
                    end
                end
            end
            
            return [n,dr_n,dtheta_n]::Array{Float64,1};
        end
        #@vectorize_2arg Float64 sys;
        idx.sys = sys;
        
        function sys!(r::Float64,theta::Float64,result::Array{Float64,1})
            result[1],result[2],result[3] = idx.sys(r,theta);
            return nothing;
        end
        idx.sys! = sys!;
        
        return idx;
    end
    
end


#General index variation from product of Chebyshev and Fourier expansions
type ChebFourIdx <: RefractiveIndex
    # ***************** TO DO: Use sparse coefficients matrix? ****************
    #coeff_1 is a 1D array with the ith entry the coefficient for the T_{i-1}(r) term.
    #Condition for smoothness at origin: coefficients for even orders of T_i(r)*cos(j*theta) or T_i(r)*sin(j*theta) must be 0. So we define:
    #coeff_c is a 2D array with the (i,j)th entry the coefficient for the T_{2*i-1}(r)*cos(j*theta) term
    #coeff_s is a 2D array with the (i,j)th entry the coefficient for the T_{2*i-1}(r)*sin(j*theta) term
    name::String;
    coeff_1::Array{Float64,1};
    coeff_c::Array{Float64,2};
    coeff_s::Array{Float64,2};
    nfunc::Function;
    sys::Function;
    sys!::Function;
    
    #Construct by fieldname
    ChebFourIdx(;coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array{Float64,0,0}) = GeneralIdx(coeff_1,coeff_c,coeff_s);
    #Construct by position
    function ChebFourIdx(coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array(Float64,(0,0)))
        idx = new("chebfour",coeff_1,coeff_c,coeff_s);
        
        #Index value
        function nfunc(r::Float64,theta::Float64)
            #Computing Chebyshev terms using recurrence relation
            n::Float64 = 0.0;
            
            #Begin Tval at T_0
            Tval_2::Float64 = 2.0*sq(r)-1; Tval_1::Float64 = r; Tval::Float64=1.0;
            
            for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
                #Compute Chebyshev value T_{i-1}(r) using recurrence relation
                Tval = 2.0*r*Tval_1 - Tval_2;
                Tval_2 = Tval_1; Tval_1 = Tval;
                
                if i <= length(idx.coeff_1)
                    n += idx.coeff_1[i]*Tval;
                end
                if mod(i,2) == 0
                    for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                        if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                            n += idx.coeff_c[i/2,j]*Tval*cos(j*theta);
                        end
                        if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                            n += idx.coeff_s[i/2,j]*Tval*sin(j*theta);
                        end
                    end
                end
            end
            return n;
        end
        @vectorize_2arg Float64 nfunc
        idx.nfunc = nfunc;
        
        #Index value and derivatives
        function sys(r::Float64,theta::Float64)
            n::Float64 = 0.0;
            dr_n::Float64 = 0.0;
            dtheta_n::Float64 = 0.0;
            
            #Begin Tval and Uval at T_0 and U_{-1} respectively
            Tval_2::Float64 = 2.0*sq(r)-1; Tval_1::Float64 = r; Tval::Float64 = 1.0;
            Uval_2::Float64 = -2.0*r; Uval_1::Float64 = -1.0; Uval::Float64 = 0.0;
            
            for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
                #Compute Chebyshev value T_{i-1}(r) using recurrence relation
                Tval = 2.0*r*Tval_1 - Tval_2;
                Tval_2 = Tval_1; Tval_1 = Tval;
                #Compute (second kind) Chebyshev value U_{i-1}(r) using recurrence relation, to compute d(T_i)/dx = i*U_{i-1}
                Uval = 2.0*r*Uval_1 - Uval_2;
                Uval_2 = Uval_1; Uval_1 = Uval;
                
                if i <= length(idx.coeff_1)
                    n += idx.coeff_1[i]*Tval;
                    dr_n += idx.coeff_1[i]*(i-1)*Uval; #Uval is one term behind
                end
                if mod(i,2) == 0
                    for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                        if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                            n += idx.coeff_c[i/2,j]*Tval*cos(j*theta);
                            dr_n += idx.coeff_c[i/2,j]*(i-1)*Uval*cos(j*theta);
                            dtheta_n -= idx.coeff_c[i/2,j]*Tval*j*sin(j*theta);
                        end
                        if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                            n += idx.coeff_s[i/2,j]*Tval*sin(j*theta);
                            dr_n += idx.coeff_s[i/2,j]*(i-1)*Uval*sin(j*theta);
                            dtheta_n += idx.coeff_s[i/2,j]*Tval*j*cos(j*theta);
                        end
                    end
                end
            end
            
            return [n,dr_n,dtheta_n]::Array{Float64,1};
        end
        #@vectorize_2arg Float64 sys
        idx.sys = sys;
        
        function sys!(r::Float64,theta::Float64,result::Array{Float64,1})
            result[1],result[2],result[3] = idx.sys(r,theta);
            return nothing;
        end
        idx.sys! = sys!;
        
        return idx;
    end
end
