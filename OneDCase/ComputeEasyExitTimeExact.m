function tau = ComputeEasyExitTimeExact(x,Bounds)
% Test: Valid only for f = 0 and g = 2

tau = -x^2/2 + (Bounds(1)+Bounds(2))/2 * x - Bounds(1)*Bounds(2)/2;

end