clear
clc
close all

% Set up of time span and space interval
dx = 0.1;
h = 0.001 * dx^2;
x = -1:dx:1;
T = 0.5;
t = 0:h:T;

% Boundary Conditions
N = length(x);
M = length(t);
u = zeros(N,M);
BoundCond = 1;
r = h / (dx^2);

% f and sigma
f = 0.1 * (32 * x.^3 - 16 * x + 1)';
sigma = 1;


if BoundCond == 0
    u = [[1;zeros(N-2,1);1],[ones(1,M-1);zeros(N-2,M-1);ones(1,M-1)]];
    
    for j = 2 : M
        u(2:end-1,j) = (1 - r*sigma^2) * u(2:end-1,j-1) + (f(2:end-1)*dx*r + 0.5*sigma^2*r) .* u(1:end-2,j-1) + (-f(2:end-1)*dx*r + 0.5*sigma^2*r) .* u(3:end,j-1);
    end
    
    figure
    plot(x,u(:,end),'o-')
elseif BoundCond == 1
    u = [[1;zeros(N-1,1)],[ones(1,M-1);zeros(N-1,M-1)]];
    for j = 2 : M
        u(2:end-1,j) = (1 - r*sigma^2) * u(2:end-1,j-1) + (f(2:end-1)*dx*r + 0.5*sigma^2*r) .* u(1:end-2,j-1) + (-f(2:end-1)*dx*r + 0.5*sigma^2*r) .* u(3:end,j-1);
        u(end,j) = u(end-1,j);
    end
    
    figure
    plot(x,u(:,end),'o-')
end
