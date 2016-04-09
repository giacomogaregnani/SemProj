function p = TwoSidedExitProbability(x0,xh,sigma,h,Bounds)

range = Bounds(2) - Bounds(1);
p = 0;

for k = -1 : 1
    p = p + exp(-2 * k * range * (k * range + xk - x0) / (sigma^2 * h)) - ...
        exp(-2 * (k * range + x0 - Bounds(2)) * (k * range + xh - Bounds(2)) / (sigma^2 * h));
end

end