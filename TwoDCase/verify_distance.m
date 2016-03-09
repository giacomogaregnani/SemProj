function check = verify_distance(X,R,M)

ind_up = find(X(2:2:end)-ones(M,1) < R);
ind_mone = find(abs(X+ones(size(X))) < R);