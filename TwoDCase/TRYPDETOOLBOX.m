model = createpde();
R1 = [3,4,-1,1,1,-1,1,1,-1,-1]';
g = decsg(R1);
geometryFromEdges(model,g);
pdegplot(g)
MSH = generateMesh(model)
Coeff = specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);
applyBoundaryCondition(model,'edge',2:4,'r',0)
applyBoundaryCondition(model,'edge',2:4,'g',0)
results = solvepde(model);
u = results.NodalSolution;
pdeplot(model,'xydata',u,'zdata',u)
PointResult = interpolateSolution(results,0,0);