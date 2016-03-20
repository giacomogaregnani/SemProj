% Find the mean exit time in the Darcy case
clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem 
Time = [0,30];
sigma = 0.01;

V = @(x,y) zeros(2,1) * x * y;
dV = @(x,y) zeros(2,1) * x * y;
f = @(x,y) -dV(x,y);
g = @(x,y) sigma * eye(2);
X0 = [0;0];
Bounds = [-1,1;-1,1];
BoundCond = 0; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
N = 2.^[3:5];
M = 1e3;

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;
phiNaive = tauNaive;
phiBernoulli = tauNaive;
tNaive = tauNaive;
tBernoulli = tauNaive;

% Solve Darcy
p = SolveDarcy(1,1);

for i = 1:length(N)
    % Compute the exit time expectation
    [tauNaive(i),phiNaive(i),tNaive(i)] = ComputeExitTimeNaiveDarcy(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,p);
    [tauBernoulli(i),phiBernoulli(i),tBernoulli(i)] = ComputeExitTimeBernoulliDarcy(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,p);
    length(N) - i
end
