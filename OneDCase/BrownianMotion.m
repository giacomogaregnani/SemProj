function W = BrownianMotion(Time,N,M)

h = (Time(2)-Time(1))/N;
W = zeros(M,N+1);
W(:,2:N+1) = cumsum(sqrt(h) * randn(M,N),2);

end
