% Define the problem 


Time = [0,100];
V = @(x) -0 * (x.^2/2 + x);
dV = @(x) -1 * (x + 1);
g = @(x) 2;
X0 = 0.5;
Bounds = [-1,1];
BoundCond = [0,0];
N = 100;
M = 10;

x = linspace(-1,1,50);

for i = 1:50
    tau_ex(i) = ComputeExitTimeExact(x(i),V,dV,g,Bounds,BoundCond);
end

plot(x,tau_ex)