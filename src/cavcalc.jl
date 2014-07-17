#require("cavbnds.jl");
require("cavidxs.jl");


#==============================================================================
#Analysis of results


#Unclusteredness heuristic measurement (lower for periodic modes)
function unclusteredness(bounces,scaling=1e4)
    #Heuristic function that characterizes the degree of unclustering in a set of points,
    #from 0 (localized on discrete points) to 1 (all over the place)
    
    if length(bounces) < 40; return (1.0, 1); end #don't evaluate short simulations
    #initiate unclusteredness score
    uclst::Float64 = 2.0;
    order::Int64 = 20;
    #check through possible orders (number of bounces in a period)
    for m = 2:20 #(for a maximum of modes with 20 bounces)
        iuclst = Array(Float64,m);
        for i = 1:m
            #Try to invent a notion of distance in modular arithmetic
            pseudomean::Float64 = angle(mean(exp(im*bounces[i:m:end,1])));
            u::Array{Float64,1} = (mod(bounces[i:m:end,1]-pseudomean+pi,2*pi)-pi)/pi;
            x::Float64 = var(u) + var(sin(bounces[i:m:end,2]));
            iuclst[i] = (1.0 - 1.0/(1+scaling*x/size(bounces,1)));
            #iuclst[i] = (1-exp(-10*x)) #Alternative function
        end
        muclst = mean(iuclst);
        if muclst < uclst
            #Record smallest
            uclst = muclst;
            order = m;
        end
    end
    return (uclst, order);
end


#Unclusteredness heuristic measurement (lower for periodic modes)
function unclusteredness2(bounces,scaling=1e4)
    #Heuristic function that characterizes the degree of unclustering in a set of points,
    #from 0 (localized on discrete points) to 1 (all over the place)
    
    if length(bounces) < 40; return (1.0, 1); end #don't evaluate short simulations
    #initiate unclusteredness score
    uclst::Float64 = 2.0;
    order::Int64 = 20;
    #check through possible orders (number of bounces in a period)
    for m = 2:20 #(for a maximum of modes with 20 bounces)
        iuclst = Array(Float64,m);
        for i = 1:m
            #Try to invent a notion of distance in modular arithmetic
            pseudomean::Float64 = angle(mean(exp(im*bounces[i:m:end,1])));
            u::Array{Float64,1} = (mod(bounces[i:m:end,1]-pseudomean+pi,2*pi)-pi)/pi;
            x::Float64 = 0.5*(var(u)+var(sin(bounces[i:m:end,2])))*m/size(bounces,1);
            iuclst[i] = (1.0 - 1.0/(1+scaling*x));
            #iuclst[i] = (1-exp(-10*x)) #Alternative function
        end
        muclst = mean(iuclst);
        if muclst < uclst
            #Record smallest
            uclst = muclst;
            order = m;
        end
    end
    return (uclst, order);
end


#Path integrator
function pathint(raypath::Array{Float64,2},idx::RefractiveIndex)
    #Calculates the integrals \int{n(s)ds} and \int{ds} over a path, which is summed 
    #to obtain the action and the path length.
    
    nval::Array{Float64,1} = nfunc(idx,raypath[:,1],raypath[:,2]);
    nval_mean::Array{Float64,1} = 0.5*(nval[2:end]+nval[1:end-1]);
    dsval::Array{Float64,1} = hypot(raypath[2:end,1].*cos(raypath[2:end,2])
                                    -raypath[1:end-1,1].*cos(raypath[1:end-1,2])
                                    ,raypath[2:end,1].*sin(raypath[2:end,2])
                                    -raypath[1:end-1,1].*sin(raypath[1:end-1,2]));
    
    pathlength::Float64 = sum(dsval);
    action::Float64 = sum(nval_mean.*dsval);
    return (pathlength,action);
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
    
    #The part of omega_real proportional to m
    dwdm::Float64 = 2*pi/action;
    #Contributions from omega_real and omega_imag independent of m
    omega0_s::Complex128 = -(sum(angle(refl_s))/cycles + im*0.5*(1 - prod(abs2(refl_s))^(1/cycles)))/action;
    omega0_p::Complex128 = -(sum(angle(refl_p))/cycles + im*0.5*(1 - prod(abs2(refl_p))^(1/cycles)))/action;
    
    return (dwdm, omega0_s, omega0_p);
end


#Compute steady state farfield emission profile
function emissionfunc(pathlengths::Array{Float64,1},phi_t1::Array{Float64,1},phi_t2::Array{Float64,1},refl_s::Array{Complex128,1},refl_p::Array{Complex128,1})
    #This assumes the cavity is in a steady state. The total output for each 
    #direction and polarization is normalized to simulate a fixed power input into 
    #each mode.
    
    assert(length(pathlengths) == length(phi_t1) == length(phi_t2) == length(refl_s) == length(refl_p));
    bouncenum::Int64 = length(pathlengths);
    
    #Transmittivities
    Trans_s::Array{Float64,1} = 1.0 - abs2(refl_s);
    Trans_p::Array{Float64,1} = 1.0 - abs2(refl_p);
    
    #Assuming uniform gain, steady state requires roundtrip gain to match loss
    Gs::Float64 = -log(prod(1.0-Trans_s))/sum(pathlengths);
    Gp::Float64 = -log(prod(1.0-Trans_p))/sum(pathlengths);
    emission_s = Array(Float64,2*bouncenum);
    emission_p = Array(Float64,2*bouncenum);
    
    #Forward direction
    power_s::Float64 = 1.0; power_p::Float64 = 1.0;
    for i = 1:bouncenum
        emission_s[i] = Trans_s[i]*power_s;
        emission_p[i] = Trans_p[i]*power_p;
        power_s *= (1.0-Trans_s[i])*exp(Gs*pathlengths[i]);
        power_p *= (1.0-Trans_p[i])*exp(Gp*pathlengths[i]);
    end
    #Normalize total output for each polarization
    emission_s[1:bouncenum] /= sum(emission_s[1:bouncenum]);
    emission_p[1:bouncenum] /= sum(emission_p[1:bouncenum]);
    
    #Reversed direction
    power_s = 1.0; power_p = 1.0;
    for i = 1:bouncenum
        power_s *= exp(Gs*pathlengths[bouncenum+1-i]);
        power_p *= exp(Gp*pathlengths[bouncenum+1-i]);
        emission_s[2*bouncenum+1-i] = Trans_s[bouncenum+1-i]*power_s;
        emission_p[2*bouncenum+1-i] = Trans_p[bouncenum+1-i]*power_p;
        power_s *= (1.0-Trans_s[bouncenum+1-i]);
        power_p *= (1.0-Trans_p[bouncenum+1-i]);
    end
    #Normalize total output for each polarization
    emission_s[bouncenum+1:end] /= sum(emission_s[bouncenum+1:end]);
    emission_p[bouncenum+1:end] /= sum(emission_p[bouncenum+1:end]);
    
    return ([phi_t1,phi_t2], emission_s, emission_p);
end


#Compute mode information
function modeeval(order::Int64, pathlengths::Array{Float64,1}, actions::Array{Float64,1}, bounces::Array{Float64,2}, bnd::Boundary, idx::RefractiveIndex)
    #Computes the frequencies and farfield pattern
    
    assert(length(pathlengths) == length(actions) == (size(bounces,1)+1));
    
    #Find integer number of cycles to analyze
    cycles::Int64 = div(size(bounces,1),order);
    bouncenum::Int64 = cycles*order;
    if bouncenum == 0
        error("Simulation time too short: For accuracy, number of bounces should be a few times larger than order of mode.");
    end
    
    #chi is angle of incidence, xi is angle of transmission
    r::Array{Float64,1}, normang::Array{Float64,1} = rsys(bnd,bounces[1:bouncenum,1]);
    n::Array{Float64,1} = nfunc(idx,r,bounces[1:bouncenum,1]); #refractive index just within edge of boundary
    coschi::Array{Float64,1} = cos(bounces[1:bouncenum,2]);
    sinxi::Array{Float64,1} = n.*sin(bounces[1:bouncenum,2]);
    
    #Fresnel reflection amplitude coefficients (for senkrecht & parallel polarizations)
    cosxi::Array{Complex128,1} = sqrt(complex128(1.0-sq(sinxi)));
    refl_s::Array{Complex128,1} = (n.*coschi-cosxi)./(n.*coschi+cosxi);
    refl_p::Array{Complex128,1} = (coschi-n.*cosxi)./(n.*cosxi+coschi);
    
    #Calculate frequency
    dwdm::Float64, omega0_s::Complex128, omega0_p::Complex128 = omegafunc(actions[2:bouncenum+1],refl_s,refl_p,cycles);
    
    #Calculate farfield emission
    notir::Array{Bool,1} = (abs(sinxi) .<= 1.0);
    phi_t1::Array{Float64,1} = nans(Float64,bouncenum);
    phi_t2::Array{Float64,1} = nans(Float64,bouncenum);
    phi_t1[notir] = normang[notir] + asin(sinxi[notir]);
    phi_t2[notir] = normang[notir] - asin(sinxi[notir]);
    phi_t::Array{Float64,1}, emission_s::Array{Float64,1}, emission_p::Array{Float64,1} = emissionfunc(pathlengths[2:bouncenum+1],phi_t1,phi_t2,refl_s,refl_p);
    
    return (dwdm, omega0_s, omega0_p, phi_t, emission_s, emission_p);
end


#Estimate mode location
function getmodebounces(order::Int64,bounces::Array{Float64,2})
    modebounces::Array{Float64,2} = Array(Float64,order,2);
    for i = 1:order
        pseudomean::Float64 = angle(mean(exp(im*bounces[i:order:end,1])));
        modebounces[i,1] = pseudomean + mean(mod(bounces[i:order:end,1]-pseudomean+pi,2*pi)-pi);
        modebounces[i,2] = mean(sin(bounces[i:order:end,2]));
    end
    return modebounces;
end


#Mode distance function
function modedist(mb1::Array{Float64,2},mb2::Array{Float64,2})
    assert(size(mb1,2) == size(mb2,2) == 2);
    mindist::Float64 = 2.0;
    if size(mb1,1) == size(mb2,1)
        order::Float64 = size(mb1,1);
        for shift = 0:order-1
            dist1::Array{Float64,1} = mb1[mod([1:order]+shift,order)+1,1]-mb2[:,1];
            dist2::Array{Float64,1} = sin(mb1[mod([1:order]+shift,order)+1,2])-sin(mb2[:,2]);
            dist::Float64 = mean(hypot(mod(abs(dist1),2*pi)/pi,dist2));
            if dist < mindist; mindist = dist; end
        end
    end
    return mindist;
end

