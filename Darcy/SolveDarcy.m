function results = SolveDarcy(sigmaA,pInlet)

% Find the solution of Darcy problem
% u = A grad P
% div(u) = 0
% with BCs to be given.

% Generate the random field
LMax = 5;
nu = 0.3;
A = realizationRF(LMax,1,nu,sigmaA,1);
NGridA = sqrt(length(A));
A = reshape(A,NGridA,NGridA);
dxA = 2 / (NGridA - 1);
[X,Y] = meshgrid(-1:dxA:1,-1:dxA:1);

% Function random field (linear interpolation)
AFunc = @(region,state) interp2(X,Y,A,region.x,region.y);

% Create PDE model
model = createpde();
R1 = [3,4,-1,1,1,-1,1,1,-1,-1]';
g = decsg(R1);
geometryFromEdges(model,g);
MSH = generateMesh(model,'Hmax',0.03);

% Specify coefficients
Coeff = specifyCoefficients(model,'m',0,'d',0,'c',1,'a',AFunc,'f',0);
applyBoundaryCondition(model,'edge',2:4,'g',0);
applyBoundaryCondition(model,'edge',3,'r',pInlet);
applyBoundaryCondition(model,'edge',1,'r',0);

results = solvepde(model);

u = results.NodalSolution;
pdeplot(model,'xydata',u,'zdata',u,'colorbar','off','colormap','default')
xlabel('x')
ylabel('y')
zlabel('p')
