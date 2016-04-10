function p = TwoSidedExitProbability(x0,xh,sigma,h,Bounds)

Sum = Bounds(2) + Bounds(1);
X = xh + x0;

if X > Sum
    p = exp(-2/(sigma^2*h) * (Bounds(2) - xh) * (Bounds(2) - x0));
else
    p = exp(-2/(sigma^2*h) * (xh - Bounds(1)) * (x0 - Bounds(1)));
end   

