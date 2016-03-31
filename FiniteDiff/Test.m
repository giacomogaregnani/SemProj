clear
clc
close all

% Set up of time span and space interval
dx = 0.1;
h = 0.01;
x = -1:dx:1;
T = 5;
t = 0:h:T;

% Boundary Conditions
N = length(x) - 2;
M = length(t);
u = zeros(N + 2,M);
BoundCond = 1;

% f and sigma
f =  @(x) -0.1 * (32 * x.^3 - 16 * x + 1)';
sigma = 1;

if BoundCond == 0
    LDiag = -(-f(x(3:end-1)) * h / (2*dx) + sigma^2 * h / (2 * dx^2));
    UDiag = -(f(x(2:end-2)) * h / (2*dx) + sigma^2 * h / (2 * dx^2));
    Diag = ones(N,1) + h / (dx^2) * sigma^2;
    A = diag(Diag) + diag(LDiag,-1) + diag(UDiag,1);
    A = sparse(A);
    BoundLeft = -f(x(2)) * h / (2 * dx) + sigma^2 * h / (2 * dx^2);
    BoundRight = f(x(end-1)) * h / (2 * dx) + sigma^2 * h / (2 * dx^2);
    
    u = [ones(1,M);zeros(N,M);ones(1,M)];
    
    for j = 2 : M
        u(2:end-1,j) = A \ [u(2,j-1) + BoundLeft; u(3:end-2,j-1); u(end-1,j-1) + BoundRight];
    end
    
elseif BoundCond == 1
    LDiag = -(-f(x(3:end)) * h / (2*dx) + sigma^2 * h / (2 * dx^2));
    UDiag = -(f(x(2:end-1)) * h / (2*dx) + sigma^2 * h / (2 * dx^2));
    Diag = ones(N+1,1) + h / (dx^2) * sigma^2;
    A = diag(Diag) + diag(LDiag,-1) + diag(UDiag,1);
    A(end,end-1) = 1;
    A(end,end) = -1;
    A = sparse(A);
    BoundLeft = -f(x(2)) * h / (2 * dx) + sigma^2 * h / (2 * dx^2);
    
    u = [ones(1,M);zeros(N+1,M)];

    for j = 2 : M
        u(2:end,j) = A \ [u(2,j-1) + BoundLeft; u(3:end-1,j-1); 0];
    end
    
end

figure
plot(x,u(:,end),'o-')




