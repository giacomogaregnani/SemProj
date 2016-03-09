function ExpTau = ComputeExitTimeBernoulliTwoD(X0,f,g,Bounds,BoundCond,W,Time)
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

if X0(1) >= Bounds(1,2) || X0(1) <= Bounds(1,1) || X0(2) >= Bounds(2,2) || X0(2) <= Bounds(2,1)
    ExpTau = 0;
    return
end

Sigma = det(g(1));
[TwoM,N] = size(W);
M = TwoM/2;
h = (Time(2)-Time(1))/(N-1);
tau = Time(2) * ones(M,1);

for j = 1:M
    w = W(2*j-1:2*j,:);
    xOld = X0;
    for i = 2:N
        xNew = EMOneStep(xOld,f,g,w(:,i)-w(:,i-1),h);
        if xNew(1) >= Bounds(1,2) || xNew(1) <= Bounds(1,1) || xNew(2) >= Bounds(2,2) || xNew(2) <= Bounds(2,1)
            tau(j) = h*(i-1);
            break
        else
            p = ComputeExitProbability(xOld,xNew,Sigma,h);
            u = rand(1,1);
            if  isempty(find(p > u,1)) == 0
                tau(j) = h*i;
                break
            end
        end
        xOld = xNew;
    end
end

ExpTau = mean(tau);

end