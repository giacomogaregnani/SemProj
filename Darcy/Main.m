% Find the mean exit time in the Darcy case
clear
close all
clc

load('ResultsRef.mat')

% Define the problem
Time = [0,10];
sigma = 0.3;
g = @(x,y) sigma * eye(2);
X0 = [-0.8; 0];
Bounds = [-1,1;-1,1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.

% Define the number of points for the grids for A
LMax = 8;

% Compute reference Random Field
% nu = 0.5;
% LC = 0.05;
% sigmaA = 1;
% A = realizationRF(LMax(end),LC,nu,sigmaA,1);
% NGridA = sqrt(length(A));
% A = reshape(A, NGridA, NGridA);

% Darcy parameters
pInlet = 1;
plotfields = 'False';

% Define the size of the grid for "preprocessing" step
deltaU = 2 .^ [-1 : -1 : -5];
J = length(deltaU);
N = Time(2) ./ deltaU;
M = 5e4;

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauCEM = zeros(1,length(N));
phiCEM = tauCEM;
tCEM = tauCEM;

for j = 1 : J
    
    % Solve Darcy
    [Ux,Uy] = SolveDarcyFF(A,pInlet,plotfields,deltaU(j));
    
    % Compute the exit time expectation
    [tauCEM(j),phiCEM(j),tCEM(j)] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(j):end),Time,Ux,Uy,deltaU(j));
    display([num2str(J - j), ' iterations remaining'])
    
    clear Ux Uy
end

% RefTau = tauCEM(end, end);
errBernoulliTau = abs(tauCEM(1 : end -1) - RefTau);

% Plot convergence
figure
loglog(deltaU(1 : end -1), errBernoulliTau,'-o')
hold on
loglog(deltaU(1 : end -1), deltaU(1 : end -1),'k--')
xlabel('h = \Delta_u')
legend('err', 'h = \Delta_u', 'Location', 'NW')
