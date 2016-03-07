% MAIN
clear
close all

% Define the problem 
Time = [0,10];
V = @(x) 0.01 * (8 * x.^4 - 8 * x.^2 + x + 2);
dV = @(x)  0.01 * (32 * x.^3 - 16 * x + 1);
f = @(x) -1 * dV(x);
g = @(x) 0.5;
X0 = 0;
Bounds = [-1,1];
BoundCond = [0,0];
N = 2.^[5:10];
M = 1000;

% Compute the BM
W = BrownianMotion(Time,N(end),M);

% Initialize
tauNaive = zeros(1,length(N));

for i = 1:length(N)
    % Compute the exit time expectation
    tauNaive(i) = ComputeExitTimeNaive(X0,f,g,Bounds,BoundCond,W(:,1:N(end)/N(i):end),Time);
    length(N) - i
end

figure
TrajectoriesForPlots(X0,f,g,Bounds,0,W,Time);

tauEx = ComputeExitTimeExact(X0,V,g,Bounds,BoundCond);
errNaive = abs(tauNaive - tauEx);

figure
loglog((Time(2)-Time(1))./N,errNaive,'ko-')
hold on
loglog((Time(2)-Time(1))./N,1./sqrt(N),'r--')

OrdersNaive = log2(errNaive(1:end-1)./errNaive(2:end));

figure
x_samp = Bounds(1):0.001:Bounds(2);
plot(x_samp,V(x_samp))