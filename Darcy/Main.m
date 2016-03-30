% Find the mean exit time in the Darcy case
clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem 
Time = [0,30];
sigma = 0.3;
g = @(x,y) sigma * eye(2);
X0 = [-0.8;-0.8];
Bounds = [-1,1;-1,1];
BoundCond = 1; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
N = 2.^[3:8];
M = 1e4;

NRef = N(end)*16;

% Compute the BM
W = BrownianMotion2D(Time,NRef,M);

% Initialize
tauNaive = zeros(1,length(N));
tauBernoulli = tauNaive;
phiNaive = tauNaive;
phiBernoulli = tauNaive;
tNaive = tauNaive;
tBernoulli = tauNaive;

% Solve Darcy
sigmaA = 1;
nu = 0.5;
LC = 0.05;
LMax = 6;
pInlet = 1;
plotfields = 'True';
[Ux,Uy,delta,tFEM] = SolveDarcy(sigmaA,LC,nu,LMax,pInlet,plotfields);

PlotVelocityField(Ux,Uy,delta)

return
for i = 1:length(N)
    % Compute the exit time expectation
    [tauNaive(i),phiNaive(i),tNaive(i)] = ComputeExitTimeNaiveDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,delta);
    [tauBernoulli(i),phiBernoulli(i),tBernoulli(i)] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time,Ux,Uy,delta);
    ComputeExitTimeNaiveDarcyPlot(X0,g,Bounds,BoundCond,W(1:30,1:N(end)/N(i):end),Time,Ux,Uy,delta);
    length(N) - i
end

save('resultsEaster')

[tauRef,phiRef,tRef] = ComputeExitTimeBernoulliDarcy(X0,g,Bounds,BoundCond,W,Time,Ux,Uy,delta);
errTauCEM = abs(tauBernoulli - tauRef);
errTauDEM = abs(tauNaive - tauRef);

save('resultsEaster')

figure
h = Time(2)./N;
loglog(h,errTauCEM,'b-*')
hold on
loglog(h,errTauDEM,'r-o')
loglog(h,h.^(0.5),'k--')
loglog(h,h,'k-')
grid on
legend('err_{CEM}','err_{DEM}','h^{0.5}','h','Location','Best')

figure
loglog(errTauCEM,tBernoulli,'b-*')
hold on
loglog(errTauDEM,tNaive,'r-o')
set(gca,'XDir','Reverse')
legend('CEM','DEM')
grid on
xlabel('error')
ylabel('work')

figure
loglog(N,tBernoulli,'b-*')
hold on
loglog(N,tNaive,'r-o')
legend('CEM','DEM')
grid on
xlabel('N')
ylabel('time')
