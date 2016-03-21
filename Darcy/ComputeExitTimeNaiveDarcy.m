function [ExpTau,ExpPhi,t] = ComputeExitTimeNaiveDarcy(X0,f,g,Bounds,BoundCond,W,Time,P,A)
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

tic
sigma = det(g(0,0));

nRand = size(A,1);
xA = linspace(-1,1,nRand);
yA = xA;
[XXA,YYA] = meshgrid(xA,yA);

if BoundCond == 0
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1) || X0(2) >= Bounds(2,2) || X0(2) <= Bounds(2,1)
        ExpTau = 0;
        return
    end
    
    [TwoM,N] = size(W);
    M = TwoM/2;
    h = (Time(2)-Time(1))/(N-1);
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    
    for j = 1:M
        w = W(2*j-1:2*j,:);
        x = X0;
        for i = 2:N
            % Evaluate the velocity field in the previous point by
            % interpolation of the results on pressure
            [ux,uy] = evaluateGradient(P,x(1),x(2));
            % Interpolate in the same point A to get the value of the
            % random variable defining the material property
            APunct = interp2(XXA,YYA,A,x(1),x(2));
            ux = -APunct * ux;
            uy = -APunct * uy;
            x = EMOneStepDarcy(x,[ux;uy],sigma,w(:,i)-w(:,i-1),h);
            if x(1) >= Bounds(1,2) || x(1) <= Bounds(1,1) || x(2) >= Bounds(2,2) || x(2) <= Bounds(2,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            end
        end
    end
    
else
    if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1)
        ExpTau = 0;
        return
    end
    
    [TwoM,N] = size(W);
    M = TwoM/2;
    h = (Time(2)-Time(1))/(N-1);
    tau = Time(2) * ones(M,1);
    phi = zeros(M,1);
    
    for j = 1:M
        w = W(2*j-1:2*j,:);
        x = X0;
        for i = 2:N
            % Evaluate the velocity field in the previous point by
            % interpolation of the results on pressure
            [ux,uy] = evaluateGradient(P,x(1),x(2));
            % Interpolate in the same point A to get the value of the
            % random variable defining the material property
            APunct = interp2(XXA,YYA,A,x(1),x(2));
            ux = -APunct * ux;
            uy = -APunct * uy;
            x = EMOneStepDarcy(x,[ux;uy],sigma,w(:,i)-w(:,i-1),h);
            if x(1) >= Bounds(1,2) || x(1) <= Bounds(1,1)
                tau(j) = h*(i-1);
                phi(j) = 1;
                break
            elseif x(2) < Bounds(2,1)
                x(2) = 2*Bounds(2,1) - x(2);
            elseif x(2) > Bounds(2,2)
                x(2) = 2*Bounds(2,2) - x(2);
            end
        end
    end
end
ExpTau = mean(tau);
ExpPhi = mean(phi);

t = toc;
end