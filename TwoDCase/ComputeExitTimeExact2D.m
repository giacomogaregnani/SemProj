function tauExact = ComputeExitTimeExact2D(N,Bounds,BoundCond,sigma,X0)

x = linspace(Bounds(1,1),Bounds(1,2),N+1);
y = linspace(Bounds(2,1),Bounds(2,2),N+1);
dx = (Bounds(2,2) - Bounds(2,1))/N;

if BoundCond == 0
    A = 1 / dx^2 * gallery('poisson',N-1);
    f = 2 / (sigma^2) * ones((N-1)^2,1);
    
    
    % Compute tau in all the space and then reshape it and impose BCs.
    tau = A \ f;
    tau = reshape(tau,N-1,N-1);
    tauFull = zeros(N+1,N+1);
    tauFull(2:end-1,2:end-1) = tau;
    
    figure
    surf(x,y,tauFull,'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    zlabel('\tau')
    
    I = find(x == X0(1));
    J = find(y == X0(2));
    
    tauExact = tauFull(I,J);
else
    % Initialize
    b = zeros((N-1)*(N+1),1);
    
    % Up and Down Boundary
    for i = 1:N-1
        A((i-1)*(N+1)+1,(i-1)*(N+1)+1:(i-1)*(N+1)+2) = h^2 * ones(1,2);
        b((i-1)*(N+1)+1) = 1;
        A(i*(N+1),i*(N+1)-1:i*(N+1)) = h^2 * ones(1,2);
        b(i*(N+1)) = 1;
    end
    
        


end