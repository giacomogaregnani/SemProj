% Find the mean exit time in the Darcy case
clear
close all
clc

% Define the problem
Time = [0,10];
sigma = 1;
g = @(x, y) sigma * eye(2);
f = @(x, y) 1 * [x^2; y];
X0 = [0; 0];
Bounds = [-1, 1; -1, 1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.

% Define the size of the grid for "preprocessing" step
deltaU = 2 .^ [-1 : -1 : -5];

J = length(deltaU);
N = 2^10;
M = 1e4;
I = length(N);

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauBernoulliC = zeros(1, length(N));
tBernoulliC = tauBernoulliC;
tauBernoulliL = zeros(1, length(N));
tBernoulliL = tauBernoulliL;

for j = 1 : J
    
    [UxConst, UyConst] = constantInterpolation(f, deltaU(j));
    [tauBernoulliC(j), ~ ,tBernoulliC(j)] = ComputeExitTimeBernoulliConst(X0,g,Bounds,BoundCond,W,Time,UxConst,UyConst,deltaU(j));
    clear UxConst UyConst
    
    [UxLin, UyLin] = linearInterpolation(f, deltaU(j));
    [tauBernoulliL(j), ~ ,tBernoulliL(j)] = ComputeExitTimeBernoulliLin(X0,g,Bounds,BoundCond,W,Time,UxLin,UyLin,deltaU(j));
    clear UxLin UyLin
    
    display([num2str(J - j), ' iterations remaining'])
    
end

[RefTau, ~, tRef] = ComputeExitTimeBernoulli2D(X0, f, g, Bounds, BoundCond, W, Time);
errBernoulliTauC = abs(tauBernoulliC - RefTau);
errBernoulliTauL = abs(tauBernoulliL - RefTau);

loglog(deltaU, errBernoulliTauC, 'o-')
hold on
loglog(deltaU, deltaU, 'k--')
grid on
xlabel('\epsilon')
legend('error', '\epsilon','Location','NW')

figure
loglog(deltaU, errBernoulliTauC, 'o-')
hold on
loglog(deltaU, errBernoulliTauL, 'o-')
loglog(deltaU, deltaU, 'k--')
loglog(deltaU, deltaU.^2, 'k')
grid on
xlabel('\epsilon')
legend('errC', 'errL', '\epsilon', '\epsilon^2')

