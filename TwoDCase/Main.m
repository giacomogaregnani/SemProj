% verify the convergence of discrete EM method for killed diffusion
% problems. 

clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem 
Time = [0,1];
sigma = 1;
V = @(x,y) zeros(2,1) * x * y;
dV = @(x,y) zeros(2,1) * x * y;
f = @(x,y) -dV(x,y);
g = @(x,y) sigma * eye(2);
X0 = [0;0];
Bounds = [-1,1;-1,1];
BoundCond = 0; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
N = 2.^[0:5];
M = 1e5;

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;
phiNaive = tauNaive;
phiBernoulli = tauNaive;

for i = 1:length(N)
    % Compute the exit time expectation
    [tauNaive(i),phiNaive(i)] = ComputeExitTimeNaive2D(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    [tauBernoulli(i),phiBernoulli(i)] = ComputeExitTimeBernoulli2D(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    length(N) - i
end

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact2D(Bounds,BoundCond,sigma,X0);
phiEx = ComputeExitProbExact2D(Bounds,BoundCond,sigma,X0,Time);
errNaivetau = abs(tauNaive - tauEx);
errBernoullitau = abs(tauBernoulli - tauEx);
errNaivephi = abs(phiNaive - phiEx);
errBernoulliphi = abs(phiBernoulli - phiEx);

% Plot the error for orders analysis on Tau
h = (Time(2)-Time(1))./N;
IndForPlots = ceil(length(N)/2);
figure
loglog(h,errNaivetau,'ro-')
hold on
loglog(h,errBernoullitau,'b*-')
loglog(h,sqrt(h)*(errNaivetau(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errBernoullitau(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('err_h^d','err_h^c','h^{0.5}','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h')

% Plot the error for orders analysis on Phi
h = (Time(2)-Time(1))./N;
IndForPlots = ceil(length(N)/2);
figure
loglog(h,errNaivephi,'ro-')
hold on
loglog(h,errBernoulliphi,'b*-')
loglog(h,sqrt(h)*(errNaivephi(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h*(errBernoulliphi(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('err_h^{d,\Phi}','err_h^{c,\Phi}','h^{0.5}','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h')

% Compute the orders
OrdersNaiveTau = log2(errNaivetau(1:end-1)./errNaivetau(2:end));
OrdersBernoulliTau = log2(errBernoullitau(1:end-1)./errBernoullitau(2:end));
OrdersNaivePhi = log2(errNaivephi(1:end-1)./errNaivephi(2:end));
OrdersBernoulliPhi = log2(errBernoulliphi(1:end-1)./errBernoulliphi(2:end));

