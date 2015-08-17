#Test RefractiveIndex abstraction
using CavChaos

idx0 = UniformIdx();
idx1 = GaussianIdx();
idx2 = TaylorIdx();
idx3 = ChebFourIdx();

#RefractiveIndex test function ================================================
function testidxfunc(idx::RefractiveIndex,rval::Array{Float64,1},thetaval::Array{Float64,1})
    println("========== Testing $(summary(idx)) ==========");
    
    #Setup test quantities
    r = Array(Float64,length(rval),length(thetaval));
    theta = Array(Float64,length(rval),length(thetaval));
    for j=1:length(thetaval)
        for i=1:length(rval)
            r[i,j] = rval[i];
            theta[i,j] = thetaval[j];
        end
    end
    
    #Scalar functions
    tic();
    n0 = similar(theta);
    n1 = similar(theta); dr_n1 = similar(theta); dtheta_n1 = similar(theta);
    for i=1:length(theta)
        n0[i] = nfunc(idx,r[i],theta[i]);
        n1[i],dr_n1[i],dtheta_n1[i] = nderiv(idx,r[i],theta[i]);
    end
    println("Scalar functions compute in $(toq())s.");
    
    #Array functions
    tic();
    n0A = nfunc(idx,r,theta);
    n1A,dr_n1A,dtheta_n1A = nderiv(idx,r,theta);
    println("Array functions compute in $(toq())s.");
    
    #Test agreement
    print("Testing function agreement... ");
    if !(n0 == n1 == n0A == n1A)
        println("n agreement failed!");
        println("Returning [r theta n0 n1 n0A n1A].");
        where = find(!(n0 .== n1 .== n0A .== n1A));
        return [r[where] theta[where] n0[where] n1[where] n0A[where] n1A[where]];
    end
    if !(dr_n1 == dr_n1A)
        println("dr_n agreement failed!");
        println("Returning [r theta dr_n1 dr_n1A].");
        where = find(!(dr_n1 .== dr_n1A));
        return [r[where] theta[where] dr_n1[where] dr_n1A[where]];
    end
    if !(dtheta_n1 == dtheta_n1A)
        println("dtheta_n agreement failed!");
        println("Returning [r theta dtheta_n1 dtheta_n1A].");
        where = find(!(dtheta_n1 .== dtheta_n1A));
        return [r[where] theta[where] dtheta_n1[where] dtheta_n1A[where]];
    end
    println("passed.");
    
    #Estimate derivatives
    dr = rval[1]; dtheta = thetaval[2];
    dr_n0 = (n0[2:end,:]-n0[1:end-1,:])/dr;
    dtheta_n0 = (n0[:,2:end]-n0[:,1:end-1])/dtheta;
    dr_n1m = 0.5*(dr_n1[1:end-1,:]+dr_n1[2:end,:]);
    dtheta_n1m = 0.5*(dtheta_n1[:,1:end-1]+dtheta_n1[:,2:end]);
    error_r = maximum(abs(dr_n0 - dr_n1m));
    error_theta = maximum(abs(dtheta_n0 - dtheta_n1m));
    
    #Test accuracy
    tol = 1e-3;
    print("Testing accuracy of derivatives... ");
    if error_r > tol
        println("dr_n error is $(string(error_r))!");
        println("Returning [r theta dr_n0 dr_n1m].");
        where = find(abs(dr_n0 - dr_n1m) .> tol);
        return [r[1:end-1,:][where] theta[1:end-1,:][where] dr_n0[where] dr_n1m[where]];
    end
    if error_theta > tol
        println("dtheta_n error is $(string(error_theta))!");
        println("Returning [r theta dtheta_n0 dtheta_n1m].");
        where = find(abs(dtheta_n0 - dtheta_n1m) .> tol);
        return [r[:,1:end-1][where] theta[:,1:end-1][where] dtheta_n0[where] dtheta_n1m[where]];
    end
    println("passed.");
    
    return "passed";
end


#Run tests ====================================================================
function testidxs()
    r = linspace(0,1.5,101)[2:end];
    theta = linspace(0,2*pi,200);

    println("############### REFRACTIVE INDEX TEST ###############");

    testinfo0 = testidxfunc(idx0,r,theta);
    testinfo1 = testidxfunc(idx1,r,theta);
    testinfo2 = testidxfunc(idx2,r,theta);
    testinfo3 = testidxfunc(idx3,r,theta);

    print("############### ");
    if testinfo0 == testinfo1 == testinfo2 == testinfo3 == "passed"
        print("PASSED ALL");
    else
        print("FAILED");
    end
    println(" ###############\n");
    
    return (testinfo0,testinfo1,testinfo2,testinfo3);
end


testidxs()