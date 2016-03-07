function ComputeExitTimeExact(X0,V,dV,g,Bounds,BoundCond)

l = Bounds(1);
r = Bounds(2);

psi = @(x) 2*((V(x) - V(l))./(g(x).^2)) .* (x-l);

FunFirstTermInside = @(z) exp(psi(z)) ./ (g(z).^2);
FirstTermInside = @(y,z) exp(-psi(y)) .* FunFirstTermInside(z);
LimSupInt = @(y) y;
FirstTerm = -2 * integral2(FirstTermInside,l,X0,l,LimSupInt);


end