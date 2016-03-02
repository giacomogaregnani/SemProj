function X = square_naive(X_0,x_bar,h,T,M,A,sigma,W)
% X = square_naive(X_0,h,T,M,A,sigma,W). 
% Vectorized version
% Finds the trajectory of a particle in the [-1,1]x[-1,1] square under a
% constant velocity field v = -Ax with the Euler-Maruyama method. 
% Input initial position X_0, timestep h, final time T, number of trajectories M
% the matrix A defining the velocity field, and stochastic factor sigma

N = ceil(T/h + 1);

X = zeros(2*M,N);
X(:,1) = repmat(X_0,M,1);
A = kron(eye(M),A);

f = @(x) -A*(x-repmat(x_bar,M,1));
g = @(x) sigma;

for i = 1:N-1
    X(:,i+1) = EM_one_step(X(:,i),f,g,W(:,i+1)-W(:,i),h);
%     index = find(abs(X(:,i+1)) > 1);
%     index_pair = index(find(mod(index,2) == 0));
%     index_odd = index(find(mod(index,2) == 1));
%     index_add_1 = index_pair - ones(size(index_pair));
%     index_add_2 = index_odd + ones(size(index_odd));
%     index = unique([index;index_add_1;index_add_2]);
end

X = X(:,1:i+1);

end