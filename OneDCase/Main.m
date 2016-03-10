% MAIN
clear
close all
clc

% Define the problem 
Time = [0,10];
V = @(x) 0.01 * (8 * x.^4 - 8 * x.^2 + x + 2);
dV = @(x)  0.01 * (32 * x.^3 - 16 * x + 1);
f = @(x) -dV(x);
g = @(x) 1;
X0 = -0.5;
Bounds = [-1,1];
BoundCond = [0,0];
N = 2.^[4:10];
M = 100;

% figure
% plot(Bounds(1):0.001:Bounds(2),V(Bounds(1):0.001:Bounds(2)))

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
IndForPlots = ceil(length(N)/2);
figure
loglog(h,errNaive,'ro-')
hold on
loglog(h,errBernoulli,'b*-')
loglog(h,sqrt(h)*(errNaive(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errBernoulli(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('err_h^d','err_h^c','h^{0.5}','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h')

% Compute the orders
OrdersNaive = log2(errNaive(1:end-1)./errNaive(2:end));
OrdersBernoulli = log2(errBernoulli(1:end-1)./errBernoulli(2:end));

% Profiles of tau vs starting point 
TauProfiles(V,dV,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(5):end),Time,10)
