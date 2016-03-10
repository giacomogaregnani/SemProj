% verify the convergence of discrete EM method for killed diffusion
% problems. 

clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem 
Time = [0,3];
sigma = 1;
V = @(x,y) zeros(2,1) * x * y;
dV = @(x,y) zeros(2,1) * x * y;
f = @(x,y) -dV(x,y);
g = @(x,y) sigma * eye(2);
X0 = [0;0];
Bounds = [-1,1;-1,1];
BoundCond = [0,0];
N = 2.^[3:12];
M = 10000;

% Compute the BM
W = BrownianMotion2D(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;

for i = 1:length(N)
    % Compute the exit time expectation
    tauNaive(i) = ComputeExitTimeNaive2D(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    tauBernoulli(i) = ComputeExitTimeBernoulli2D(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    length(N) - i
end

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact2D(30,Bounds,sigma,X0);
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

