require("cavbnds.jl");
import GSL;

#==============================================================================
#List of refractive index distributions


#Each refractive index distribution idx has at position (r,theta) a refractive index n = nfunc(idx,r,theta), with derivatives dr_n and dtheta_n for the radial and angular components respectively. nsys(idx,r,theta) returns an array containing n, dr_n and dtheta_n, while nsys!(idx,r,theta,A) stores the result in array A.

abstract RefractiveIndex;

function nsys!(idx::RefractiveIndex,r::Float64,theta::Float64,result::Array{Float64,1})
    result[1],result[2],result[3] = nsys(idx,r,theta);
    return nothing;
end


#Uniform ----------------------------------------------------------------------
type UniformIdx <: RefractiveIndex
    #1 argument: the refractive index
    index_base::Float64
    
    #Construct by fieldname
    UniformIdx(;index_base::Float64=2.0) = new(index_base);
    #Construct by abbreviation
    UniformIdx(n0::Float64=2.0) = new(n0);
    #Construct by position
    UniformIdx(n0::Float64=2.0) = new(n0);
end
        
#Index value
nfunc(idx::UniformIdx,r::Float64,theta::Float64) = idx.index_base::Float64;
function nfunc{N}(idx::UniformIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta));
    n::Array{Float64,N} = fill(idx.index_base,size(r));
    return n::Array{Float64,N};
end

#Index value and derivatives
nsys(idx::UniformIdx,r::Float64,theta::Float64) = (idx.index_base,0.0,0.0);
function nsys{N}(idx::UniformIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta));
    return (fill(idx.index_base,size(r)), zeros(Float64,size(r)), zeros(Float64,size(r)));
end


#Gaussian ---------------------------------------------------------------
type GaussianIdx <: RefractiveIndex
    #4 arguments: the base (underlying uniform) refractive index, and the remaining
    #index_variation, spread and center fields are arrays where each row 
    #corresponds to the center positions, amplitude and amount of spread (stddev) of 
    #a Gaussian variation contribution. Each row of the spread has 3 components: 2 
    #for the standard deviations in 2 perpendicular directions and the last for the 
    #rotation angle. The overall index is the sum of the base index and the gaussian 
    #variations.
    index_base::Float64;
    center::Array{Float64,2};
    index_variation::Array{Float64,1};
    spread::Array{Float64,2};
    
    #Construct by fieldname
    GaussianIdx(;index_base::Float64=2.0,center::Array{Float64,2}=[0.0 0.0],index_variation::Array{Float64,1}=[-0.2],spread::Array{Float64,2}=[0.5 0.5 0.0]) = GaussianIdx(index_base,center,index_variation,spread);
    #Construct by abbreviation
    GaussianIdx(;n0::Float64=2.0,c::Array{Float64,2}=[0.0 0.0],dn::Array{Float64,1}=[-0.2],sig::Array{Float64,2}=[0.5 0.5 0.0]) = GaussianIdx(n0,c,dn,sig);
    #Construct by position
    function GaussianIdx(n0::Float64=2.0,c::Array{Float64,2}=[0.0 0.0],dn::Array{Float64,1}=[-0.2],sig::Array{Float64,2}=[0.5 0.5 0.0])
        if !(length(dn) == size(sig,1) == size(c,1))
            error("Index variation, spread and center arrays must be of the same number of rows.");
        end
        return new(n0,c,dn,sig);
    end
end
        
#Index value
function nfunc(idx::GaussianIdx,r::Float64,theta::Float64)
    #assert(length(idx.index_variation) == size(idx.center,1) == size(idx.spread,1));
    x::Float64 = r*cos(theta);
    y::Float64 = r*sin(theta);
    n::Float64 = idx.index_base;
    
    for i = 1:length(idx.index_variation)
        dxr::Float64 = (x-idx.center[i,1])*cos(idx.spread[i,3]) - (y-idx.center[i,2])*sin(idx.spread[i,3]);
        dyr::Float64 = (y-idx.center[i,2])*cos(idx.spread[i,3]) + (x-idx.center[i,1])*sin(idx.spread[i,3]);
        n += idx.index_variation[i]*exp(-0.5*(sq(dxr/idx.spread[i,1]) + sq(dyr/idx.spread[i,2])));
    end
    return n::Float64;
end

function nfunc{N}(idx::GaussianIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(length(idx.index_variation) == size(idx.center,1) == size(idx.spread,1));
    assert(size(r)==size(theta));
    x::Array{Float64,N} = r.*cos(theta);
    y::Array{Float64,N} = r.*sin(theta);
    n::Array{Float64,N} = fill(idx.index_base,size(r));
    
    for i = 1:length(idx.index_variation)
        dxr::Array{Float64,N} = (x-idx.center[i,1])*cos(idx.spread[i,3]) - (y-idx.center[i,2])*sin(idx.spread[i,3]);
        dyr::Array{Float64,N} = (y-idx.center[i,2])*cos(idx.spread[i,3]) + (x-idx.center[i,1])*sin(idx.spread[i,3]);
        n += idx.index_variation[i]*exp(-0.5*(sq(dxr/idx.spread[i,1]) + sq(dyr/idx.spread[i,2])));
    end
    return n::Array{Float64,N};
end

#Index value and derivatives
function nsys(idx::GaussianIdx,r::Float64,theta::Float64)
    #assert(length(idx.index_variation) == size(idx.center,1) == size(idx.spread,1));
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
    return (n,dr_n,dtheta_n);
end

function nsys{N}(idx::GaussianIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(length(idx.index_variation) == size(idx.center,1) == size(idx.spread,1));
    assert(size(r)==size(theta));
    x::Array{Float64,N} = r.*cos(theta);
    y::Array{Float64,N} = r.*sin(theta);
    n::Array{Float64,N} = fill(idx.index_base,size(r));
    dr_n::Float64 = zeros(Float64,size(r));
    dtheta_n::Float64 = zeros(Float64,size(r));
    
    for i = 1:length(idx.index_variation)
        dxr::Array{Float64,N} = (x-idx.center[i,1])*cos(idx.spread[i,3]) - (y-idx.center[i,2])*sin(idx.spread[i,3]);
        dyr::Array{Float64,N} = (y-idx.center[i,2])*cos(idx.spread[i,3]) + (x-idx.center[i,1])*sin(idx.spread[i,3]);
        gaussians::Array{Float64,N} = idx.index_variation[i]*exp(-0.5*(sq(dxr/idx.spread[i,1]) + sq(dyr/idx.spread[i,2])));
        
        n += gaussians;
        dr_n -= (cos(theta+idx.spread[i,3]).*dxr/sq(idx.spread[i,1]) + sin(theta+idx.spread[i,3]).*dyr/sq(idx.spread[i,2])).*gaussians;
        dtheta_n += ((y*cos(idx.spread[i,3])+x*sin(idx.spread[i,3])).*dxr/sq(idx.spread[i,1]) + (-x*cos(idx.spread[i,3])+y*sin(idx.spread[i,3])).*dyr/sq(idx.spread[i,2])).*gaussians;
    end
    return (n,dr_n,dtheta_n);
end


#Taylor -----------------------------------------------------------------------
type TaylorIdx <: RefractiveIndex
    # ***************** TO DO: Use sparse coefficients matrix? ****************
    #General index variation from Taylor expansion
    #1 argument: coeff is a matrix where the (i,j)th element is the coefficient of 
    #the x^(i-1)*y^(j-1) term in the Taylor expansion about the origin of the index 
    #distribution
    coeff::Array{Float64,2};
    
    #Construct by fieldname
    TaylorIdx(;coeff::Array{Float64,2}=[2.0 0.0;0.0 0.4]) = new(coeff);
    #Construct by abbreviation
    TaylorIdx(;cf::Array{Float64,2}=[2.0 0.0;0.0 0.4]) = new(cf);
    #Construct by position
    TaylorIdx(coeff::Array{Float64,2}=[2.0 0.0;0.0 0.4]) = new(coeff);
end
        
#Index value
function nfunc(idx::TaylorIdx,r::Float64,theta::Float64)
    x::Float64 = r*cos(theta);
    y::Float64 = r*sin(theta);
    n::Float64 = 0.0;
    for i = 1:size(idx.coeff,1)
        for j = 1:size(idx.coeff,2)
            n += idx.coeff[i,j]*GSL.pow_int(y,i-1)*GSL.pow_int(x,j-1);
        end
    end
    return n::Float64;
end

function nfunc{N}(idx::TaylorIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta));
    x::Array{Float64,N} = r.*cos(theta);
    y::Array{Float64,N} = r.*sin(theta);
    n::Array{Float64,N} = zeros(Float64,size(r));
    for i = 1:size(idx.coeff,1)
        for j = 1:size(idx.coeff,2)
            for k = 1:length(r)
                n[k] += idx.coeff[i,j]*GSL.pow_int(y[k],i-1)*GSL.pow_int(x[k],j-1);
            end
        end
    end
    return n::Array{Float64,N};
end

#Index value and derivatives
function nsys(idx::TaylorIdx,r::Float64,theta::Float64)
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
    return (n,dr_n,dtheta_n);
end

function nsys{N}(idx::TaylorIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta));
    x::Array{Float64,N} = r.*cos(theta);
    y::Array{Float64,N} = r.*sin(theta);
    n::Array{Float64,N} = zeros(Float64,size(r));
    dr_n::Array{Float64,N} = zeros(Float64,size(r));
    dtheta_n::Array{Float64,N} = zeros(Float64,size(r));
    
    for i = 1:size(idx.coeff,1)
        for j = 1:size(idx.coeff,2)
            for k = 1:length(r)
                n[k] += idx.coeff[i,j]*GSL.pow_int(y[k],i-1)*GSL.pow_int(x[k],j-1);
                if i >= 2
                    dr_n[k] += idx.coeff[i,j]*(i-1)*GSL.pow_int(y[k],i-2)*GSL.pow_int(x[k],j-1)*sin(theta[k]);
                    dtheta_n[k] += idx.coeff[i,j]*(i-1)*GSL.pow_int(y[k],i-2)*GSL.pow_int(x[k],j);
                end
                if j >= 2
                    dr_n[k] += idx.coeff[i,j]*GSL.pow_int(y[k],i-1)*(j-1)*GSL.pow_int(x[k],j-2)*cos(theta[k]);
                    dtheta_n[k] -= idx.coeff[i,j]*GSL.pow_int(y[k],i)*(j-1)*GSL.pow_int(x[k],j-2);
                end
            end
        end
    end
    return (n,dr_n,dtheta_n);
end


#Chebyshev-Fourier ------------------------------------------------------------
type ChebFourIdx <: RefractiveIndex
    # ***************** TO DO: Use sparse coefficients matrix? ****************
    #General index variation from product of Chebyshev and Fourier expansions
    #3 arguments:
    #coeff_1 is a 1D array with the ith entry the coefficient for the T_{i-1}(r) term.
    #Condition for smoothness at origin: coefficients for even orders of T_i(r)*cos(j*theta) or T_i(r)*sin(j*theta) must be 0. So we define:
    #coeff_c is a 2D array with the (i,j)th entry the coefficient for the T_{2*i-1}(r)*cos(j*theta) term
    #coeff_s is a 2D array with the (i,j)th entry the coefficient for the T_{2*i-1}(r)*sin(j*theta) term
    coeff_1::Array{Float64,1};
    coeff_c::Array{Float64,2};
    coeff_s::Array{Float64,2};
    
    #Construct by fieldname
    ChebFourIdx(;coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array{Float64,0,0}) = new(coeff_1,coeff_c,coeff_s);
    #Construct by abbreviation
    ChebFourIdx(;c1::Array{Float64,1}=[2.0],cc::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],cs::Array{Float64,2}=Array{Float64,0,0}) = new(c1,cc,cs);
    #Construct by position
    ChebFourIdx(coeff_1::Array{Float64,1}=[2.0],coeff_c::Array{Float64,2}=[0.1 -0.1; -0.1 0.1],coeff_s::Array{Float64,2}=Array(Float64,(0,0))) = new(coeff_1,coeff_c,coeff_s);
end
        
#Index value
function nfunc(idx::ChebFourIdx,r::Float64,theta::Float64)
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
    return n::Float64;
end

function nfunc{N}(idx::ChebFourIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta));
    #Computing Chebyshev terms using recurrence relation
    n::Array{Float64,N} = zeros(Float64,size(r));
    
    #Begin Tval at T_0
    Tval_2::Array{Float64,N} = 2.0*sq(r)-1.0; Tval_1::Array{Float64,N} = r; Tval::Array{Float64,N}=ones(Float64,size(r));
    
    for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
        #Compute Chebyshev value T_{i-1}(r) using recurrence relation
        Tval = 2.0*r.*Tval_1 - Tval_2;
        Tval_2 = Tval_1; Tval_1 = Tval;
        
        if i <= length(idx.coeff_1)
            n += idx.coeff_1[i]*Tval;
        end
        if mod(i,2) == 0
            for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                    n += idx.coeff_c[i/2,j]*Tval.*cos(j*theta);
                end
                if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                    n += idx.coeff_s[i/2,j]*Tval.*sin(j*theta);
                end
            end
        end
    end
    return n::Array{Float64,N};
end
        
#Index value and derivatives
function nsys(idx::ChebFourIdx,r::Float64,theta::Float64)
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
    return (n,dr_n,dtheta_n);
end

function nsys{N}(idx::ChebFourIdx,r::Array{Float64,N},theta::Array{Float64,N})
    assert(size(r)==size(theta));
    n::Array{Float64,N} = zeros(Float64,size(r));
    dr_n::Array{Float64,N} = zeros(Float64,size(r));
    dtheta_n::Array{Float64,N} = zeros(Float64,size(r));
    
    #Begin Tval and Uval at T_0 and U_{-1} respectively
    Tval_2::Array{Float64,N} = 2.0*sq(r)-1.0; Tval_1::Array{Float64,N} = r; Tval::Array{Float64,N} = ones(Float64,size(r));
    Uval_2::Array{Float64,N} = -2.0*r; Uval_1::Array{Float64,N} = fill(-1.0,size(r)); Uval::Array{Float64,1} = zeros(Float64,size(r));
    
    for i = 1:max(length(idx.coeff_1),2*size(idx.coeff_c,1),2*size(idx.coeff_s,1))
        #Compute Chebyshev value T_{i-1}(r) using recurrence relation
        Tval = 2.0*r.*Tval_1 - Tval_2;
        Tval_2 = Tval_1; Tval_1 = Tval;
        #Compute (second kind) Chebyshev value U_{i-1}(r) using recurrence relation, to compute d(T_i)/dx = i*U_{i-1}
        Uval = 2.0*r.*Uval_1 - Uval_2;
        Uval_2 = Uval_1; Uval_1 = Uval;
        
        if i <= length(idx.coeff_1)
            n += idx.coeff_1[i]*Tval;
            dr_n += idx.coeff_1[i]*(i-1)*Uval; #Uval is one term behind
        end
        if mod(i,2) == 0
            for j = 1:max(size(idx.coeff_c,2),size(idx.coeff_s,2))
                if i <= 2*size(idx.coeff_c,1) && j <= size(idx.coeff_c,2)
                    n += idx.coeff_c[i/2,j]*Tval.*cos(j*theta);
                    dr_n += idx.coeff_c[i/2,j]*(i-1)*Uval.*cos(j*theta);
                    dtheta_n -= idx.coeff_c[i/2,j]*j*Tval.*sin(j*theta);
                end
                if i <= 2*size(idx.coeff_s,1) && j <= size(idx.coeff_s,2)
                    n += idx.coeff_s[i/2,j]*Tval.*sin(j*theta);
                    dr_n += idx.coeff_s[i/2,j]*(i-1)*Uval.*sin(j*theta);
                    dtheta_n += idx.coeff_s[i/2,j]*j*Tval.*cos(j*theta);
                end
            end
        end
    end
    return (n,dr_n,dtheta_n);
end

