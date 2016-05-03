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

nStep = 0;

if BoundCond == 0
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1) || X0(2) >= Bounds(2,2) || X0(2) <= Bounds(2,1)
        ExpTau = 0;
        return
    end
    
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    sigma = g(1,1);
    sigma = sigma(1,1);
    nStep = phi;
    
    for j = 1:M
        xOld = X0;
        time = 0;
        while time < Time(2)
            
            xNew = EMOneStep(xOld,f,sigma,h);
            nStep(j) = nStep(j) + 1;
            time = time + h;
            
            if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1) || xNew(2) >= Bounds(2,2) || xNew(2) <= Bounds(2,1)
                tau(j) = time;
                phi(j) = 1;
                time = Time(2);
            else
                p = ComputeExitProbability(xOld,xNew,sigma,h);
                u = rand(1,1);
                
                if  isempty(find(p > u,1)) == 0
                    tau(j) = time;
                    phi(j) = 1;
                    time = Time(2);
                end
                
            end
        end
    end
    
elseif BoundCond == 1
    
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1)
        ExpTau = 0;
        return
    end
    
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    sigma = g(1,1);
    sigma = sigma(1,1);
    nStep = phi;
    
    for j = 1:M
        xOld = X0;
        time = 0;
        while time < Time(2)
            
            if time + h > Time(2)
                h = Time(2) - time;
            end
            
            xNew = EMOneStep(xOld,f,sigma,h);
            nStep(j) = nStep(j) + 1;
            
            time = time + h;
            
            if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1)
                tau(j) = time;
                phi(j) = 1;
                time = Time(2);
                
            elseif xNew(2) < Bounds(2,1)
                xNew(2) = 2*Bounds(2,1) - xNew(2);
                
            elseif xNew(2) > Bounds(2,2)
                xNew(2) = 2*Bounds(2,2) - xNew(2);
                
            else
                
                p = ComputeExitProbability(xOld,xNew,sigma,h);
                u = rand(1,1);
                p(2) = 0;
                p(4) = 0;
                
                if  isempty(find(p > u,1)) == 0
                    tau(j) =time;
                    phi(j) = 1;
                    time = Time(2);
                end
            end
        end
    end
end

ExpTau = mean(tau);
ExpPhi = mean(phi);
nStep = mean(nStep);
end