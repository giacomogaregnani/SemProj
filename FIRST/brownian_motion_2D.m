function W = brownian_motion_2D(t_0,T,h,M)

N = ceil((T-t_0)/h + 1);
W = zeros(2*M,N);
W(:,2:N) = cumsum(sqrt(h) * randn(2*M,N-1),2);

end