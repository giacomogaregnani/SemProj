function [Pressure,ux,uy] = SolveDarcy(sigmaA,pInlet,plotfields)

% Find the solution of Darcy problem
% u = A grad P
% div(u) = 0
% with BCs to be given.

% Generate the random field
LMax = 5;
nu = 0.1;
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
applyBoundaryCondition(model,'edge',[1,3],'g',0);
applyBoundaryCondition(model,'edge',4,'r',pInlet);
applyBoundaryCondition(model,'edge',2,'r',0);

results = solvepde(model);

AInt = zeros(size(model.Mesh.Nodes,2),1);
p = model.Mesh.Nodes;

for i = 1 : length(AInt)
    AInt(i) = interp2(X,Y,A,p(1,i),p(2,i));
end

ux = -results.XGradients .* AInt;
uy = -results.YGradients .* AInt;

if strcmp(plotfields,'True') == 1
    Pressure = results.NodalSolution;
    pdeplot(model,'xydata',Pressure,'zdata',Pressure,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('p')
    
    figure
    pdeplot(model,'xydata',AInt,'zdata',AInt,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('A')
    
    figure
    pdeplot(model,'xydata',ux,'zdata',ux,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('u_x')
    
    figure
    pdeplot(model,'xydata',uy,'zdata',uy,'colorbar','off','colormap','default')
    xlabel('x')
    ylabel('y')
    zlabel('u_y')
end


