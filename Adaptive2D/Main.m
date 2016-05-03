% verify the convergence of DEM and CEM in a square domain.

clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% Define the problem
Time = [0, 3];
sigma = 1;

f = @(x,y) 0 * [x; y];
g = @(x,y) sigma * eye(2);
X0 = [0; 0];
Bounds = [-1,1;-1,1];
BoundCond = 0; % 0 for killing everywhere. 1 for two killing and two reflecting BCs.
l = 0 : 7;
h = (Time(2) - Time(1)) ./ (2.^l);
hmin = 2.^(-l) .* h;
M = 1e4;

% Initialize
tauNaiveAd = zeros(1,length(h));
phiNaiveAd = tauNaiveAd;
tNaiveAd = tauNaiveAd;
tauNaive = tauNaiveAd;
phiNaive = tauNaiveAd;
tNaive = tNaiveAd;
tauCEM = tauNaiveAd;
phiCEM = tauNaiveAd;
tCEM = tauNaiveAd;

tauEx = ComputeExitTimeExact2D(Bounds,BoundCond,sigma,X0);

for j = 1:length(h)
    % Compute the exit time expectation
    [tauNaiveAd(j), phiNaiveAd(j), tNaiveAd(j), hDEMAd, DDemAd] = DEM2DAdapt(X0, f, g, Bounds, BoundCond, M, Time, hmin(j), h(j));
    disp('Ad Completed')
    [tauNaive(j), phiNaive(j), tNaive(j)] = DEM2D(X0, f, g, Bounds, BoundCond, M, Time, hmin(j));
    disp('DEM Completed')
    [tauCEM(j), phiCEM(j), tCEM(j)] = CEM2D(X0, f, g, Bounds, BoundCond, M, Time, h(j));
    disp('CEM Completed')
    disp([num2str(length(h) - j), ' iterations remaining'])
end

% Compute the exact expectation of tau and the error
errNaivetauAd = abs(tauNaiveAd - tauEx);
errNaivetau =  abs(tauNaive - tauEx);
errCEM = abs(tauCEM - tauEx);

% Plot the error for orders analysis on Tau
IndForPlots = ceil(length(h)/2);
figure
loglog(h, errNaivetauAd,'ro-')
hold on
loglog(h, errNaivetau,'b-*')
loglog(h, errCEM, '-<')
loglog(h,h*(errNaivetauAd(IndForPlots)/h(IndForPlots)),'k')
grid on
h_legend = legend('adaptive','DEM, h_{bound}','CEM, h_{int}', 'h_{int}');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('h_{int}')
ylabel('error')

% plot number of steps
figure
loglog(h, tNaiveAd, 'r-o')
hold on
loglog(h, tNaive, 'b-*')
loglog(h, tCEM, '-<')
grid on
xlabel('h_{int}')
ylabel('Mean number of timesteps')
h_legend = legend('adaptive','DEM, h_{bound}', 'CEM, h_{int}')
set(h_legend,'Location','northeast','FontSize',13);


