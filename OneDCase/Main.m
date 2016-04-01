% Program to compute the orders of convergence of DEM and CEM methods.

% MAIN
clear
close all
clc

% Define the problem 
Time = [0,1];
% SMOOTH f
V = @(x) 0 * x.^0;
dV = V;
% V = @(x) 1 * (8 * x.^4 - 8 * x.^2 + x + 2);
% dV = @(x)  1 * (32 * x.^3 - 16 * x + 1);
% ROUGH f
% V = @(x) 0.1 * ((-1 - 2*x) .* (x < -0.5) + (4*x + 2) .* (x >= -0.5) .* (x < 0) + (2 - 2*x) .* (x >= 0) .* (x < 0.5) + (4*x - 1) .* (x >= 0.5));
% dV = @(x) 0.1 * (-2 * (x < -0.5) + 4 * (x >= -0.5) .* (x < 0) -2 * (x >= 0) .* (x < 0.5) + 4 * (x >= 0.5));
% f = @(x) -dV(x);
f = @(x) 0 * x.^0;
g = @(x) 1.5;
X0 = 0;
Bounds = [-1,1];
BoundCond = [0,1];
N = 2.^[0:12];
M = 1e4;

% figure
% plot(Bounds(1):0.001:Bounds(2),V(Bounds(1):0.001:Bounds(2)),'LineWidth',2)
% hold on
% plot(Bounds(1):0.001:Bounds(2),dV(Bounds(1):0.001:Bounds(2)),'-.r','LineWidth',2)
% plot(Bounds(1):0.001:Bounds(2),VRough(Bounds(1):0.001:Bounds(2)),'--g','LineWidth',2)
% plot(Bounds(1):0.001:Bounds(2),dVRough(Bounds(1):0.001:Bounds(2)),':k','LineWidth',2)
% legend('V smooth','dV smooth','V rough','dV rough','Location','NW')
% grid on

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
% TrajectoriesForPlots(X0,f,g,Bounds,0,W,Time);

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact(X0,V,g,Bounds,BoundCond);
phiEx = ComputeExitProbFD(X0,Time,Bounds,BoundCond,f,g(1));
errNaive = abs(tauNaive - tauEx);
errBernoulli = abs(tauBernoulli - tauEx);
errNaivePhi = abs(phiNaive - phiEx);
errBernoulliPhi = abs(phiBernoulli - phiEx);

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

figure
loglog(h,errNaivePhi,'ro-')
hold on
loglog(h,errBernoulliPhi,'b*-')
loglog(h,sqrt(h)*(errNaivePhi(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errBernoulliPhi(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('err_h^d','err_h^c','h^{0.5}','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h')

% Computational time
% figure
% loglog(errNaive,tNaive,'r-o')
% hold on
% loglog(errBernoulli,tBern,'b-*')
% h_legend = legend('time_{DEM}','time_{CEM}');
% set(h_legend,'Location','northwest','FontSize',13);
% set(gca,'XDir','Reverse')
% xlabel('error')
% ylabel('computational time')
% grid on

% Compute the orders
OrdersNaive = log2(errNaive(1:end-1)./errNaive(2:end));
OrdersBernoulli = log2(errBernoulli(1:end-1)./errBernoulli(2:end));
OrdersNaivePhi = log2(errNaivePhi(1:end-1)./errNaivePhi(2:end));
OrdersBernoulliPhi = log2(errBernoulliPhi(1:end-1)./errBernoulliPhi(2:end));


% % Profiles of tau vs starting point 
TauProfiles(V,dV,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(7):end),Time,10)
PhiProfiles(f,g,Bounds,BoundCond,W(1:1000,1:N(end)/N(7):end),Time,10);
