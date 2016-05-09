% Find the mean exit time in the Darcy case
clear
close all
clc

% Define the problem
Time = [0,10];
sigma = 1;
g = @(x,y) sigma * eye(2);
X0 = [-0.8;-0.8];
Bounds = [-1,1;-1,1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
N = 2.^12;
M = 1;

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;
phiNaive = tauNaive;
phiBernoulli = tauNaive;
tNaive = tauNaive;
tBernoulli = tauNaive;

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
deltaU = 2 ./ [2 : 16]; % 2 .^ [0 : -1 : -6, -11];

J = length(deltaU);
I = length(N);

for j = 1 : J
    
    % Solve Darcy
    [Ux,Uy] = SolveDarcyFF(A,pInlet,plotfields,deltaU(j));
    
    for i = 1 : I
        % Compute the exit time expectation
%         [tauNaive(j,i),phiNaive(j,i),tNaive(j,i)] = ComputeExitTimeNaiveDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,delta(j));
        [tauBernoulli(j,i),phiBernoulli(j,i),tBernoulli(j,i)] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,deltaU(j));
%         ComputeExitTimeNaiveDarcyPlot(X0,g,Bounds,BoundCond,W(1:30,1:N(end)/N(i):end),Time,Ux,Uy,delta);
        display([num2str(I*J-((j-1)*I + i)), ' iterations remaining'])
    end
    
    clear Ux Uy    
end

RefTau = tauBernoulli(end, end);
errBernoulliTau = abs(tauBernoulli - RefTau);
errNaiveTau = abs(tauNaive - RefTau);
ConvergencePlots(errBernoulliTau, errNaiveTau, tauBernoulli, deltaU, N, Time);