% MAIN

% Define the problem 
Time = [0,100];
f = @(x) -0.1 * (x + 0.5);
g = @(x) 0.05;
X0 = 0.5;
Bounds = [-1,1];
N = 100;
M = 10;

% Compute the BM
W = BrownianMotion(Time,N,M);

% Compute the exit time expectation
tau = ComputeExitTimeNaive(X0,f,g,Bounds,0,W,Time);

TrajectoriesForPlots(X0,f,g,Bounds,0,W,Time);