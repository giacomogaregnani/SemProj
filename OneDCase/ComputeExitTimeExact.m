function ComputeExitTimeExact(X0,f,g,Bounds,BoundCond)

l = Bounds(1);
r = Bounds(2);

psi = @(x) integral(f,l,x);

FunFirstTermInside = @(z) exp(psi(z)) ./ (g(z).^2);
FirstTermInside = @(y) exp(-psi(y)) .* (y/NInt) .* cumtrapz(FunFirstTermInside(linspace(0,y,NInt));
FirstTerm = -2 * integral(FirstTermInside,l,X0);

end