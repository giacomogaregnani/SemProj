function [ExpTau, ExpPhi, nStep] = CEM2D(X0, f, g, Bounds, BoundCond, M, Time, h)
% ExpTau = ComputeExitTimeBernoulli(X0,f,g,Bounds,BoundCond,N,M)
% Compute expected exit time with Euler-Maruyama method with Bernoulli
% implementation of the killed boundary condition.
% Input : X0 starting point; f,g functions defining the equation dX =
% f(X)dt + g(X)dW; Bounds a vector [l,u] of the lower and upper bound
% for the equations; BoundCond a vector [b1,b2] defining the boundary
% conditions (bi = 0 killing, bi = 1 reflecting); W a matrix containing M
% realisations of a one-dimensional BM on N intervals in the time-span [t0,T];
% Time the vector [t0,T]

% Bounds(1,:) for x-direction, Bounds(2,:) for y direction

if BoundCond == 0
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1) || X0(2) >= Bounds(2,2) || X0(2) <= Bounds(2,1)
        ExpTau = 0;
        return
    end
    
    Sigma = g(1,1);
    Sigma = Sigma(1,1);
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    nStep = phi;
    N = Time(2) / h;
    
    for j = 1:M
        xOld = X0;
        for i = 2:N
            xNew = EMOneStep(xOld,f,Sigma,h);
            nStep(j) = nStep(j) + 1;
            
            if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1) || xNew(2) >= Bounds(2,2) || xNew(2) <= Bounds(2,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            else
                p = ComputeExitProbability(xOld,xNew,Sigma,h);
                u = rand(1,1);
                if  isempty(find(p > u,1)) == 0
                    tau(j) = h*(i-1);
                    phi(j) = 1;
                    break
                end
            end
            xOld = xNew;
        end
    end
    
else
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1)
        ExpTau = 0;
        return
    end
    
    Sigma = g(1,1);
    Sigma = Sigma(1,1);
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    nStep = phi;
    N = Time(2) / h;
    
    for j = 1:M
        xOld = X0;
        for i = 2:N
            xNew = EMOneStep(xOld,f,Sigma,h);
            nStep(j) = nStep(j) + 1;
            
            if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1) 
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            elseif xNew(2) < Bounds(2,1)
                xNew(2) = 2*Bounds(2,1) - xNew(2);
            elseif xNew(2) > Bounds(2,2)
                xNew(2) = 2*Bounds(2,2) - xNew(2);
            else
                p = ComputeExitProbability(xOld,xNew,Sigma,h);
                p = [p(1),0,p(3),0];
                u = rand(1,1);
                if  isempty(find(p > u,1)) == 0
                    tau(j) = h*(i-1);
                    phi(j) = 1;
                    break
                end
            end
            xOld = xNew;
        end
    end
end

ExpTau = mean(tau);
ExpPhi = mean(phi);
nStep = mean(nStep);

end
