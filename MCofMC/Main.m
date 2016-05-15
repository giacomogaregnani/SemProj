% Find the mean exit time in the Darcy case
clear
close all
clc

% Domain and general parameters
Time = [0,20];
sigma = [1, 0.5];
initialCondition = [-0.8; 0];
bounds = [-1,1;-1,1];
boundCond = 1;

% Parameters of outer MC
nRealizations = 10;
deltaU = 2 .^ -3;
pInlet = 1;
plotFields = 'False';
LMax = 8;
nu = 0.5;
LC = 0.05;
sigmaA = 1;

% Parameters of inner MC
N = Time(2) ./ deltaU;
nTrajCEM = 1e4;

% Initialization
nSigma = length(sigma);
tauCEM = zeros(1, nSigma);

for j = 1 : nSigma
    for i = 1 : nRealizations
        
        % Generate the random field
        A = realizationRF(LMax(end),LC,nu,sigmaA,1);
        NGridA = sqrt(length(A));
        A = reshape(A, NGridA, NGridA);
        
        % Solve Darcy
        [Ux,Uy] = SolveDarcyFF(A,pInlet,plotFields,deltaU);
        
        % Generate a BM
        W = BrownianMotion2D(Time, N, nTrajCEM);
        
        % Compute the exit time expectation
        currTauCEM = ComputeExitTimeBernoulliDarcy(initialCondition, sigma(j), bounds, boundCond, W, Time, Ux, Uy, deltaU);
        tauCEM(j) = tauCEM(j) + currTauCEM;
        
        display([num2str(nRealizations - i), ' iterations remaining'])
    end
    tauCEM(j) = tauCEM(j) / nRealizations;
end
