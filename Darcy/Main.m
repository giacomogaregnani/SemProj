% Find the mean exit time in the Darcy case
clear
close all
clc

% Define the problem
Time = [0,100];
sigma = 0.3;
g = @(x,y) sigma * eye(2);
X0 = [-0.8;-0.8];
Bounds = [-1,1;-1,1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.

% Define the number of points for the grids for A
LMax = 8;

% Compute reference Random Field
nu = 0.5;
LC = 0.05;
sigmaA = 1;
A = realizationRF(LMax(end),LC,nu,sigmaA,1);
NGridA = sqrt(length(A));
A = reshape(A, NGridA, NGridA);

% Darcy parameters
pInlet = 1;
plotfields = 'False';

% Define the size of the grid for "preprocessing" step
deltaU = 2 .^ [-1 : -1 : -5, -12];
J = length(deltaU);
N = 2^10;
M = 1e5;
I = length(N);

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauBernoulli = zeros(1,length(N));
phiBernoulli = tauBernoulli;
tBernoulli = tauBernoulli;

for j = 1 : J
    
    % Solve Darcy
    [Ux,Uy] = SolveDarcyFF(A,pInlet,plotfields,deltaU(j));
    
    for i = 1 : I
        % Compute the exit time expectation
        [tauBernoulli(j,i),phiBernoulli(j,i),tBernoulli(j,i)] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,deltaU(j));
        display([num2str(I*J-((j-1)*I + i)), ' iterations remaining'])
    end
    
    clear Ux Uy    
end

RefTau = tauBernoulli(end, end);
errBernoulliTau = abs(tauBernoulli(1 : end -1) - RefTau);

% Plot convergence wrt deltaU
loglog(deltaU(1 : end -1), errBernoulliTau,'-o')
hold on
loglog(deltaU(1 : end -1), deltaU(1 : end -1),'k--')
xlabel('\Delta_u')
legend('err', '\Delta_u', 'Location', 'NW')

% Misc : Convergence plots.
% ConvergencePlots(errBernoulliTau, errNaiveTau, tauBernoulli, deltaU, N, Time);