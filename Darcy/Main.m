% Find the mean exit time in the Darcy case
clear
close all
clc

% Define the problem
Time = [0,10];
sigma = 0.3;
g = @(x,y) sigma * eye(2);
X0 = [-0.8;-0.8];
Bounds = [-1,1;-1,1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
N = 2.^[0 : 5];
M = 1e5;

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
LMax = [1, 4, 7];

% Compute reference Random Field
nu = 0.5;
LC = 0.05;
sigmaA = 1;
A = realizationRF(LMax(end),LC,nu,sigmaA,1);
NGridA = sqrt(length(A));
A = reshape(A,NGridA,NGridA);

% Darcy parameters
pInlet = 1;
plotfields = 'F';

J = length(LMax);
I = length(N);
delta = zeros(J,1);

for j = 1 : J
    
    % Solve Darcy
    restrX = 1 : 2 ^ (LMax(end) - LMax(j)) : size(A,1);
    [Ux,Uy,delta(j)] = SolveDarcyFF(A(restrX,restrX),pInlet,plotfields);
    
    for i = 1 : I
        % Compute the exit time expectation
%         [tauNaive(j,i),phiNaive(j,i),tNaive(j,i)] = ComputeExitTimeNaiveDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,delta(j));
        [tauBernoulli(j,i),phiBernoulli(j,i),tBernoulli(j,i)] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,delta(j));
%         ComputeExitTimeNaiveDarcyPlot(X0,g,Bounds,BoundCond,W(1:30,1:N(end)/N(i):end),Time,Ux,Uy,delta);
        display([num2str(I*J-((j-1)*I + i)), ' iterations remaining'])
    end
    
    clear Ux Uy    
end

RefTau = tauBernoulli(end, end);
errBernoulliTau = abs(tauBernoulli - RefTau);
errNaiveTau = abs(tauNaive - RefTau);
ConvergencePlots(errBernoulliTau, errNaiveTau, tauNaive, tauBernoulli, delta, N, Time);