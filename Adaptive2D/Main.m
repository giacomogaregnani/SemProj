% verify the convergence of DEM and CEM in a square domain.

clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem
Time = [0,2];
sigma = 1;

f = @(x,y) 0 * [x; y];
g = @(x,y) sigma * eye(2);
X0 = [0;0];
Bounds = [-1,1;-1,1];
BoundCond = 0; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
l = 5 : 7;
h = (Time(2) - Time(1)) ./ (2.^l);
M = 1e4;

% Initialize
tauNaiveAd = zeros(1,length(h));
phiNaiveAd = tauNaiveAd;
tNaiveAd = tauNaiveAd;
tauNaive = tauNaiveAd;
phiNaive = tauNaiveAd;
tNaive = tNaiveAd;

for j = 1:length(h)
    % Compute the exit time expectation
    [tauNaiveAd(j), phiNaiveAd(j), tNaiveAd(j), hDEMAd] = DEM2DAdapt(X0, f, g, Bounds, BoundCond, M, Time, h(j), h(1));
    [tauNaive(j), phiNaive(j), tNaive(j)] = DEM2D(X0, f, g, Bounds, BoundCond, M, Time, h(j));
    
    disp(length(h) - j)
end

% Compute the exact expectation of tau and the error
tauEx = ComputeExitTimeExact2D(Bounds,BoundCond,sigma,X0);
phiEx = ComputeExitProbExact2D(Bounds,BoundCond,sigma,X0,Time);
errNaivetauAd = abs(tauNaiveAd - tauEx);
errNaivephiAd = abs(phiNaiveAd - phiEx);
errNaivetau =  abs(tauNaive - tauEx);
errNaivephi = abs(phiNaive - phiEx);

% Plot the error for orders analysis on Tau
IndForPlots = ceil(length(h)/2);
figure
loglog(h,errNaivetauAd,'ro-')
hold on
loglog(h,sqrt(h)*(errNaivetauAd(IndForPlots)/sqrt(h(IndForPlots))),'k--')
loglog(h,h,'k')
grid on
h_legend = legend('err_h^{d,\tau}','h^{0.5}','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h')

% Plot the error for orders analysis on Phi
figure
loglog(h,errNaivephiAd,'ro-')
hold on
loglog(h,sqrt(h),'k--')
loglog(h,h,'k')
grid on
h_legend = legend('err_h^{d,\Phi}','h^{0.5}','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h')

% plot err vs time
figure
loglog(errNaivetau,tNaive,'b--*')
hold on
loglog(errNaivetauAd,tNaiveAd,'r--o')
grid on


% Compute the orders
OrdersNaiveTau = log2(errNaivetauAd(1:end-1)./errNaivetauAd(2:end));
OrdersNaivePhi = log2(errNaivephiAd(1:end-1)./errNaivephiAd(2:end));

