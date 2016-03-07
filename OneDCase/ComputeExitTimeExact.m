function tau = ComputeExitTimeExact(X0,V,g,Bounds,BoundCond)

l = Bounds(1);
r = Bounds(2);

psi = @(x) 2 * ( V(l) - V(x) ) ./ (g(x).^2); % Formula is valid only if g is constant and f = -V'.

FirstTermArg = @(x,y) exp(-psi(x)) .* exp(psi(y)) ./ (g(y).^2);
LimSupInt = @(x) x;
FirstTerm = integral2(FirstTermArg,l,X0,l,LimSupInt);

SecondTermArg = @(y) exp(-psi(y));
SecondTerm = integral(SecondTermArg,l,X0);

if BoundCond(2) == 0
    c1 = 2 * integral2(FirstTermArg,l,r,l,LimSupInt) / integral(SecondTermArg,l,r);
else
    c1Arg = @(z) exp(psi(z)) ./ (g(z).^2);
    c1 = 2 * integral(c1Arg,l,r);
end

tau = -2 * FirstTerm + c1 * SecondTerm;

end