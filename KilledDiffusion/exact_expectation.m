function [X,exp_out,exp_tau] = exact_expectation(A,sigma,X_0,x_bar,T,h,M,W)
% Compute the analytical solution of the SDE in the square for an initial
% position X_0 and an attractor x_bar

N = ceil(T/h);
X = zeros(2*M,N+1);
a_x = A(1,1);
a_y = A(2,2);

t_samp = 0:h:T;

X(:,1) = repmat(X_0,M,1);


% compute solution in x direction
t_matrix_x_det = repmat(exp(-a_x*t_samp),M,1);
t_matrix_x_sto = repmat(exp(a_x*t_samp),M,1);
deterministic_part_x = x_bar(1) * ones(M,N) + t_matrix_x_det(:,2:N+1).*((X_0(1) - x_bar(1))*ones(M,N));
X(1:2:end-1,2:N+1) = deterministic_part_x + sigma * t_matrix_x_det(:,2:end) .* cumsum((W(1:2:end-1,2:end)-W(1:2:end-1,1:end-1)) .* t_matrix_x_sto(:,2:end),2);

% compute solution in y direction
t_matrix_y_det = repmat(exp(-a_y*t_samp),M,1);
t_matrix_y_sto = repmat(exp(a_y*t_samp),M,1);
deterministic_part_y = x_bar(2) * ones(M,N) + t_matrix_y_det(:,2:N+1).*((X_0(2) - x_bar(2))*ones(M,N));
X(2:2:end,2:N+1) = deterministic_part_y + sigma * t_matrix_y_det(:,2:end) .* cumsum((W(2:2:end,2:end)-W(2:2:end,1:end-1)) .* t_matrix_y_sto(:,2:end),2);

out = zeros(M,1);
tau = out;
 
for i = 1:M
    out(i) = 1 - isempty(find(abs(X(2*i-1,:)) > 1,1)) * isempty(find(abs(X(2*i,:)) > 1,1));
    if isempty(find(abs(X(2*i-1,:)) > 1,1)) && isempty(find(abs(X(2*i,:)) > 1,1))
        tau(i) = T;
    elseif isempty(find(abs(X(2*i-1,:)) > 1,1))
        tau(i) = min(T,find(abs(X(2*i,:)) > 1,1)*h);
    elseif isempty(find(abs(X(2*i,:)) > 1,1))
        tau(i) = min(T,find(abs(X(2*i-1,:)) > 1,1)*h);
    else
        tau(i) = min(find(abs(X(2*i,:)) > 1,1)*h,find(abs(X(2*i-1,:)) > 1,1)*h);
    end
end

exp_out = sum(out)/M;
exp_tau = sum(tau)/M;
end