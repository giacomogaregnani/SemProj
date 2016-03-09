function [exp_out,exp_tau] = square_naive_memory(X_0,x_bar,h,T,M,A,sigma,W)
% X = square_bernoulli(X_0,h,T,M,A,sigma,W). 


N = ceil(T/h + 1);

X = repmat(X_0,M,1);

f = @(x) -A*(x-x_bar);
g = @(x) sigma;
tau = T * ones(M,1);
out = zeros(M,1);

for j = 1:M
    for i = 1:N-1
        X(2*j-1:2*j) = EM_one_step(X(2*j-1:2*j),f,g,W(2*j-1:2*j,i+1)-W(2*j-1:2*j,i),h);
      
        if abs(X(2*j-1)) > 1 || abs(X(2*j)) > 1
            tau(j) = h*i;
            out(j) = 1;
            break
        end
    end
end

exp_out = mean(out);
exp_tau = mean(tau);

end