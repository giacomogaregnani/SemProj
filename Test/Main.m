% Find the mean exit time in the Darcy case
clear
close all
clc

% Define the problem
Time = [0,10];
sigma = 1;
g = @(x, y) sigma * eye(2);
f = @(x, y) [abs(x); y];
X0 = [0; 0];
Bounds = [-1, 1; -1, 1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.

% Define the size of the grid for "preprocessing" step
deltaU = 2 .^ [-1 : -1 : -8];

J = length(deltaU);
N = 2^10;
M = 1e4;
I = length(N);

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauBernoulli = zeros(1,length(N));
phiBernoulli = tauBernoulli;
tBernoulli = tauBernoulli;

for j = 1 : J
    
    [Ux, Uy] = constantInterpolation(f, deltaU(j));
    
    for i = 1 : I
        % Compute the exit time expectation
        [tauBernoulli(j,i), ~ ,tBernoulli(j,i)] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,deltaU(j));
        display([num2str(I*J-((j-1)*I + i)), ' iterations remaining'])
    end
   
    clear Ux Uy
end

[RefTau, ~, tRef] = ComputeExitTimeBernoulli2D(X0, f, g, Bounds, BoundCond, W, Time);
errBernoulliTau = abs(tauBernoulli - RefTau);

loglog(deltaU(1 : end - 1), errBernoulliTau(1 : end - 1), 'o-')
hold on
loglog(deltaU(1 : end - 1), deltaU(1 : end - 1))
grid on
xlabel('\epsilon')
legend('err', '\epsilon')

