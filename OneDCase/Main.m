% MAIN
clear
close all
clc

% Define the problem 
Time = [0,10];
V = @(x) 0.01 * (8 * x.^4 - 8 * x.^2 + x + 2);
dV = @(x)  0.01 * (32 * x.^3 - 16 * x + 1);
f = @(x) -dV(x);
g = @(x) 3;
X0 = 0.5;
Bounds = [-1,1];
BoundCond = [0,0];
N = 2.^[5:13];
M = 1000;

% Compute the BM
W = BrownianMotion(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;

for i = 1:length(N)
    % Compute the exit time expectation
    tauNaive(i) = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    tauBernoulli(i) = ComputeExitTimeBernoulli(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    length(N) - i
end

% Plot the trajectories for the finest timestep (Comment for memory and CPUT save)
% figure
% TrajectoriesForPlots(X0,f,g,Bounds,0,W,Time);

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact(X0,V,g,Bounds,BoundCond);
errNaive = abs(tauNaive - tauEx);
errBernoulli = abs(tauBernoulli - tauEx);

% Plot the error for orders analysis
h = (Time(2)-Time(1))./N;
figure
loglog(h,errNaive,'ro-')
hold on
loglog(h,sqrt(h)*(errNaive(6)/sqrt(h(6))),'k--')
loglog(h,errBernoulli,'bo-')
loglog(h,h*(errBernoulli(6)/h(6)),'k--')
grid on

% Compute the orders
OrdersNaive = log2(errNaive(1:end-1)./errNaive(2:end));
OrdersBernoulli = log2(errBernoulli(1:end-1)./errBernoulli(2:end));

% figure
% x_samp = Bounds(1):0.001:Bounds(2);
% plot(x_samp,V(x_samp))