#Test Boundary abstraction
using CavChaos

bnd0 = CircularBnd();
bnd1 = EllipticalBnd();
bnd2 = PolygonalBnd(5,0.5);
bnd3 = FourierBnd();


#Boundary test function
function testbndfunc(bnd::Boundary,theta::Array{Float64,1})
    println("========== Testing $(summary(bnd)) ==========");
    
    #Scalar functions
    tic();
    r0 = similar(theta);
    r1 = similar(theta); dtheta_r1 = similar(theta);
    r2 = similar(theta); dtheta_r2 = similar(theta); dtheta2_r2 = similar(theta);
    for i=1:length(theta)
        r0[i] = rfunc(bnd,theta[i]);
        r1[i],dtheta_r1[i] = rderiv(bnd,theta[i]);
        r2[i],dtheta_r2[i],dtheta2_r2[i] = rderiv2(bnd,theta[i]);
    end
    println("Scalar functions compute in $(toq())s.");
    
    #Array functions
    tic();
    r0A = rfunc(bnd,theta);
    r1A,dtheta_r1A = rderiv(bnd,theta);
    r2A,dtheta_r2A,dtheta2_r2A = rderiv2(bnd,theta);
    println("Array functions compute in $(toq())s.");
    
    #Test agreement
    print("Testing function agreement... ");
    if !(r0 == r1 == r2 == r0A == r1A == r2A)
        println("r agreement failed!");
        println("Returning [theta r0 r1 r2 r0A r1A r2A].");
        where = find(!(r0 .== r1 .== r2 .== r0A .== r1A .== r2A));
        return [theta[where] r0[where] r1[where] r2[where] r0A[where] r1A[where] r2A[where]];
    end
    if !(dtheta_r1 == dtheta_r2 == dtheta_r1A == dtheta_r2A)
        println("dtheta_r agreement failed!");
        println("Returning [theta dtheta_r1 dtheta_r2 dtheta_r1A dtheta_r2A].");
        where = find(!(dtheta_r1 .== dtheta_r2 .== dtheta_r1A .== dtheta_r2A));
        return [theta[where] dtheta_r1[where] dtheta_r2[where] dtheta_r1A[where] dtheta_r2A[where]];
    end
    if !(dtheta2_r2 == dtheta2_r2A)
        println("dtheta2_r agreement failed.");
        println("Returning [theta dtheta2_r2 dtheta2_r2A].");
        where = find(!(dtheta2_r2 .== dtheta2_r2A));
        return [theta[where] dtheta2_r2[where] dtheta2_r2A[where]];
    end
    println("passed.");
    
    #Estimate derivatives
    dtheta = theta[2];
    dtheta_r0 = (r0[2:end]-r0[1:end-1])/dtheta; #length-1
    dtheta2_r0 = (dtheta_r0[2:end]-dtheta_r0[1:end-1])/dtheta; #length-2
    dtheta2_r1 = (dtheta_r1[2:end]-dtheta_r1[1:end-1])/dtheta; #length-1
    dtheta_r1m = 0.5*(dtheta_r1[1:end-1]+dtheta_r1[2:end]);
    dtheta2_r2m = 0.5*(dtheta2_r2[1:end-1]+dtheta2_r2[2:end]);
    if typeof(bnd) == PolygonalBnd
        #exclude corners for PolygonalBnd
        angdiv = 2*pi/bnd.sides;
        sides = (abs(mod(theta+0.5*angdiv,angdiv)-0.5*angdiv) .> dtheta)
        error_0 = maximum(abs(dtheta_r0 - dtheta_r1m)[sides[1:end-1]]);
        error2_0 = maximum(abs(dtheta2_r0 - dtheta2_r2[2:end-1])[sides[2:end-1]]);
        error2_1 = maximum(abs(dtheta2_r1 - dtheta2_r2m)[sides[1:end-1]]);
    else
        error_0 = maximum(abs(dtheta_r0 - dtheta_r1m));
        error2_0 = maximum(abs(dtheta2_r0 - dtheta2_r2[2:end-1]));
        error2_1 = maximum(abs(dtheta2_r1 - dtheta2_r2m));
    end
    
    #Test accuracy
    tol = 1e-4;
    print("Testing accuracy of derivatives... ");
    if error_0 > tol
        println("dtheta_r1 error is $(string(error_0))!");
        println("Returning [theta dtheta_r0 dtheta_r1m].");
        where = find(abs(dtheta_r0 - dtheta_r1m) .> tol);
        return [theta[where] dtheta_r0[where] dtheta_r1m[where]];
    end
    if error2_0 > tol
        println("dtheta2_r2 error is $(tring(error2_0))!");
        println("Returning [theta dtheta2_r0 dtheta2_r2].");
        where = find(abs(dtheta2_r0 - dtheta2_r2[2:end-1]) .> tol);
        return [theta[2:end-1][where] dtheta2_r0[where] dtheta2_r2[2:end-1][where]];
    end
    if error2_1 > tol
        println("dtheta2_r2 error is $(string(error2_1))!");
        println("Returning [theta dtheta2_r1 dtheta2_r2m].");
        where = find(abs(dtheta2_r1 - dtheta2_r2m) .> tol);
        return [theta[where] dtheta2_r1[where] dtheta2_r2m[where]];
    end
    println("passed.");
    
    return "passed";
end


function testbnds()
    theta = linspace(0,2*pi,1000);

    println("############### BOUNDARY TEST ###############");

    testinfo0 = testbndfunc(bnd0,theta);
    testinfo1 = testbndfunc(bnd1,theta);
    testinfo2 = testbndfunc(bnd2,theta);
    testinfo3 = testbndfunc(bnd3,theta);

    print("############### ");
    if testinfo0 == testinfo1 == testinfo2 == testinfo3 == "passed"
        print("PASSED ALL");
    else
        print("FAILED");
    end
    println(" ###############\n");
    
    return (testinfo0,testinfo1,testinfo2,testinfo3);
end


testbnds()