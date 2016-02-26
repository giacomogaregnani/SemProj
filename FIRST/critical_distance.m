function bool = critical_distance(r,tol)
% i-th boolean is 1 if the particle is too near to i-th wall.

bool = 1 - (r-tol > 0);

end