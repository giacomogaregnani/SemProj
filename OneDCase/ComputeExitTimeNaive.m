function ExpTau = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,W,Time)
% ExpTau = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,N,M)
% Compute expected exit time with Euler-Maruyama method with naive
% implementation of the killed boundary condition.
% Input : X0 starting point; f,g functions defining the equation dX =
% f(X)dt + g(X)dW; Bounds a vector [l,u] of the lower and upper bound
% for the equations; BoundCond a vector [b1,b2] defining the boundary
% conditions (bi = 0 killing, bi = 1 reflecting); W a matrix containing M 
% realisations of a one-dimensional BM on N intervals in the time-span [t0,T];
% Time the vector [t0,T]

[M,N] = size(W);
h = (Time(2)-Time(1))/N;
tau = Time(2) * ones(M,1);

x = X0;

for j = 1:M
    w = W(j,:);
    x = X0;
    for i = 2:N
        x = EMOneStep(x,f,g,w(i)-w(i-1),h);
        if x > Bounds(2) || x < Bounds(1)
            tau(j) = h*(i-1);
            break
        end
    end
end

ExpTau = mean(tau);
end