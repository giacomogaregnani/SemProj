function [exp_out,exp_tau] = square_bernoulli(X_0,x_bar,h,T,M,A,sigma,W,tol)
% X = square_bernoulli(X_0,h,T,M,A,sigma,W). 


N = ceil(T/h + 1);

X = zeros(2*M,N);
X(:,1) = repmat(X_0,M,1);

f = @(x) -A*(x-x_bar);
g = @(x) sigma;
tau = T * ones(M,1);
out = zeros(M,1);

for j = 1:M
    for i = 1:N-1
        X(2*j-1:2*j,i+1) = EM_one_step(X(2*j-1:2*j,i),f,g,W(2*j-1:2*j,i+1)-W(2*j-1:2*j,i),h);
      
        if abs(X(2*j-1,i+1)) > 1 || abs(X(2*j,i+1)) > 1
            tau(j) = h*i;
            out(j) = 1;
            break
        else
            [r_1,r_2,r_3,r_4] = compute_distance(X(2*j-1:2*j,i+1));
            bool = critical_distance([r_1,r_2,r_3,r_4],tol);
            p = compute_exit_probability(X(2*j-1:2*j,i),X(2*j-1:2*j,i+1),bool,sigma,h);
            if  isempty(find(p > 0.5,1)) == 0
                tau(j) = h*i;
                out(j) = 1;
                break
            end
        end
    end
end

exp_out = mean(out);
exp_tau = mean(tau);

end