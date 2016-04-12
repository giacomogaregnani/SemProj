% Program to compute the orders of convergence of DEM and CEM methods.

% MAIN
clear
close all
clc

% Define the problem
Time = [0,0.5];
% SMOOTH f
V = @(x) 1 * (8 * x.^4 - 8 * x.^2 + x + 2);
dV = @(x) 1 * (32 * x.^3 - 16 * x + 1);
% ROUGH f
% V = @(x) 0.1 * ((-1 - 2*x) .* (x < -0.5) + (4*x + 2) .* (x >= -0.5) .* (x < 0) + (2 - 2*x) .* (x >= 0) .* (x < 0.5) + (4*x - 1) .* (x >= 0.5));
% dV = @(x) 0.1 * (-2 * (x < -0.5) + 4 * (x >= -0.5) .* (x < 0) -2 * (x >= 0) .* (x < 0.5) + 4 * (x >= 0.5));
f = @(x) -dV(x);
g = @(x) 2;
X0 = 0;
Bounds = [-1,1];
BoundCond = [0,0];
N = 2.^[3 : 12];
M = 1e3;

for k = 1 : 1
    % Compute the BM
    W = BrownianMotion(Time,N(end),M);
    
    % Initialize
    tauNaive = zeros(1,length(N));
    tauBernoulli = tauNaive;
    tNaive = tauNaive;
    tBern = tauBernoulli;
    phiNaive = tauNaive;
    phiBernoulli = tauNaive;
    
    for i = 1:length(N)
        % Compute the exit time expectation
        [tauNaive(i),phiNaive(i),tNaive(i)] = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
        [tauBernoulli(i),phiBernoulli(i),tBern(i)] = ComputeExitTimeBernoulli(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
        length(N) - i
    end
    
    % Plot the trajectories for the finest timestep (Comment for memory and CPUT save)
    % figure
    % TrajectoriesForPlots(X0,f,g,Bounds,0,W(1:10,1:N(end)/N(3):end),Time);
    
    % Compute the exact expectation of tau and the error
    tauEx = ComputeExitTimeExact(X0,V,g,Bounds,BoundCond);
    phiEx = ComputeExitProbFD2(X0,Time,Bounds,BoundCond,f,g(1));
    
    errNaive(k,:) = abs(tauNaive - tauEx);
    errBernoulli(k,:) = abs(tauBernoulli - tauEx);
    errNaivePhi(k,:) = abs(phiNaive - phiEx);
    errBernoulliPhi(k,:) = abs(phiBernoulli - phiEx);
    
end

if k > 1
    errNaive = mean(errNaive);
    errNaivePhi = mean(errNaivePhi);
    errBernoulli = mean(errBernoulli);
    errBernoulliPhi = mean(errBernoulliPhi);
end

ErrorPlots(errBernoulliPhi,errNaivePhi,errBernoulli,errNaive,N,Time)

% Compute the orders
OrdersNaive = log2(errNaive(1:end-1)./errNaive(2:end));
OrdersBernoulli = log2(errBernoulli(1:end-1)./errBernoulli(2:end));
OrdersNaivePhi = log2(errNaivePhi(1:end-1)./errNaivePhi(2:end));
OrdersBernoulliPhi = log2(errBernoulliPhi(1:end-1)./errBernoulliPhi(2:end));

% % Profiles of tau vs starting point
% TauProfiles(V,dV,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(end-1):end),Time,10);
% PhiProfiles(f,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(end-2):end),Time,10);
