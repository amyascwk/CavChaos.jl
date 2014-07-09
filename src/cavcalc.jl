#require("cavbnds.jl");
require("cavidxs.jl");


#==============================================================================
#Analysis of results


#Path integrator
function pathint(raypath::Array{Float64,2},idx::RefractiveIndex)
    #Calculates the integrals \int{n(s)ds} and \int{ds} over a path, which is summed 
    #to obtain the action and the path length.
    
    nval::Array{Float64,1} = idx.nfunc(raypath[:,1],raypath[:,2]);
    nval_mean::Array{Float64,1} = 0.5*(nval[2:end]+nval[1:end-1]);
    
    x1val = raypath[1:end-1,1].*cos(raypath[1:end-1,2]);
    y1val = raypath[1:end-1,1].*sin(raypath[1:end-1,2]);
    x2val = raypath[2:end,1].*cos(raypath[2:end,2]);
    y2val = raypath[2:end,1].*sin(raypath[2:end,2]);
    dsval = hypot(x2val-x1val,y2val-y1val);
    
    pathlength::Float64 = sum(dsval);
    action::Float64 = sum(nval_mean.*dsval);
    return (pathlength,action);
end


#Evaluate bounces data to compute Fresnel coefficients and farfield emission profile
function bounceseval(bounces::Array{Float64,2},bnd::Boundary,idx::RefractiveIndex)
    #Calculates the Fresnel reflection amplitude coefficients (for both senkrecht and parallel polarizations)
    #Calculates the farfield directions and Transmittivities of emitted rays
    
    #chi is angle of incidence, xi is angle of transmission
    theta::Array{Float64,1} = bounces[:,1];
    chi::Array{Float64,1} = bounces[:,2];
    r::Array{Float64,1}, normang::Array{Float64,1} = bnd.sys(theta);
    n::Array{Float64,1} = idx.nfunc(r,theta); #refractive index just within edge of boundary
    coschi::Array{Float64,1} = cos(chi);
    sinxi::Array{Float64,1} = n.*sin(chi);
    phi_t1::Array{Float64,1} = nans(Float64,size(bounces,1));
    phi_t2::Array{Float64,1} = nans(Float64,size(bounces,1));
    phi_t1[sinxi.<=1.0] = normang[sinxi.<=1.0] + asin(sinxi[sinxi.<=1.0]);
    phi_t2[sinxi.<=1.0] = normang[sinxi.<=1.0] - asin(sinxi[sinxi.<=1.0]);
    cosxi::Array{Complex128,1} = sqrt(complex128(1.0-sq(sinxi)));
    
    #Reflection coefficients
    refl_s::Array{Complex128,1} = (n.*coschi-cosxi)./(n.*coschi+cosxi);
    refl_p::Array{Complex128,1} = (coschi-n.*cosxi)./(n.*cosxi+coschi);
    
    #Transmittivities
    Trans_s::Array{Float64,1} = 1.0 - abs2(refl_s);
    Trans_p::Array{Float64,1} = 1.0 - abs2(refl_p);
    
    return (refl_s, refl_p, phi_t1, phi_t2, Trans_s, Trans_p);
end


#Compute allowed complex frequencies
function omegafunc(actions::Array{Float64,1},refl_s::Array{Complex128,1},refl_p::Array{Complex128,1},cycles::Int64=1)
    #Returns parameters to calculate the mth allowed complex angular frequency 
    #omega_m/c for a periodic mode, where
    #omega_m/c = omega_real - im*omega_imag = dwdm*m + omega0
    #in the dimensionless length units used in cavbnds.jl and cavidxs.jl
    
    #Total phase around the loop is \Phi = \oint{k(s)ds} + \sum{\arg refl} = omega_real*action + \sum{\arg refl} = 2*pi*m
    #Fractional loss per cycle is dUdt/U*T = dUdt/U*\oint{ds/(c/n(s))} = (-2*omega_imag*c)*(action/c) = 1-\prod{|refl|^2}
    
    assert(length(actions) == length(refl_s) == length(refl_p));
    action::Float64 = sum(actions)/cycles;
    
    dwdm::Float64 = 2*pi/action;
    omega0_s::Complex128 = -(sum(angle(refl_s))/cycles + im*0.5*(1 - prod(abs2(refl_s))^(1/cycles)))/action;
    omega0_p::Complex128 = -(sum(angle(refl_p))/cycles + im*0.5*(1 - prod(abs2(refl_p))^(1/cycles)))/action;
    
    return (dwdm, omega0_s, omega0_p);
end


#Compute steady state farfield emission profile
function emissionfunc(pathlengths::Array{Float64,1},phi_t1::Array{Float64,1},phi_t2::Array{Float64,1},Trans_s::Array{Float64,1},Trans_p::Array{Float64,1})
    #This assumes the cavity is in a steady state.
    
    assert(length(pathlengths) == length(phi_t1) == length(phi_t2) == length(Trans_s) == length(Trans_p));
    bouncenum::Int64 = length(pathlengths);
    
    #Assuming uniform gain
    Gs::Float64 = -log(prod(1.0-Trans_s))/sum(pathlengths);
    Gp::Float64 = -log(prod(1.0-Trans_p))/sum(pathlengths);
    output_s = Array(Float64,2*bouncenum);
    output_p = Array(Float64,2*bouncenum);
    
    #Forward direction
    power_s::Float64 = 1.0;
    power_p::Float64 = 1.0;
    for i = 1:bouncenum
        output_s[i] = Trans_s[i]*power_s;
        output_p[i] = Trans_p[i]*power_p;
        power_s *= (1.0-Trans_s[i])*exp(Gs*pathlengths[i]);
        power_p *= (1.0-Trans_p[i])*exp(Gp*pathlengths[i]);
    end
    
    #Reversed direction
    #************* Issue: How to normalize output emission, in order to compare both 
    #rays travelling in both directions, and compare different cavities? ***********
    power_s = 1.0;
    power_p = 1.0;
    for i = 1:bouncenum
        power_s *= exp(Gs*pathlengths[bouncenum+1-i]);
        power_p *= exp(Gp*pathlengths[bouncenum+1-i]);
        output_s[2*bouncenum+1-i] = Trans_s[bouncenum+1-i]*power_s;
        output_p[2*bouncenum+1-i] = Trans_p[bouncenum+1-i]*power_p;
        power_s *= (1.0-Trans_s[bouncenum+1-i]);
        power_p *= (1.0-Trans_p[bouncenum+1-i]);
    end
    
    return ([phi_t1,phi_t2], output_s, output_p);
end


#Compute mode information
function modeeval(order::Int64, pathlengths::Array{Float64,1}, actions::Array{Float64,1}, bounces::Array{Float64,2}, bnd::Boundary, idx::RefractiveIndex)
    #Computes the frequencies and farfield pattern
    
    assert(length(pathlengths) == length(actions) == (size(bounces,1)+1));
    cycles::Int64 = div(size(bounces,1),order);
    bouncenum::Int64 = cycles*order;
    
    if bouncenum == 0
        error("Number of bounces should be a few times larger than order of mode.");
    else
        #Important quantities
        refl_s::Array{Complex128,1}, refl_p::Array{Complex128,1}, phi_t1::Array{Float64,1}, phi_t2::Array{Float64,1}, Trans_s::Array{Float64,1}, Trans_p::Array{Float64,1} = bounceseval(bounces[1:bouncenum,:],bnd,idx);
        
        #Frequency (let's set m = 100 for now)
        dwdm::Float64, omega0_s::Complex128, omega0_p::Complex128 = omegafunc(actions[2:bouncenum+1],refl_s,refl_p,cycles);
        
        #Farfield emission
        phi_t::Array{Float64,1}, output_s::Array{Float64,1}, output_p::Array{Float64,1} = emissionfunc(pathlengths[2:bouncenum+1],phi_t1,phi_t2,Trans_s,Trans_p);
    end
    
    return (dwdm, omega0_s, omega0_p, phi_t, output_s, output_p);
end