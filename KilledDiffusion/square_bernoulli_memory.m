function [exp_out,exp_tau] = square_bernoulli_memory(X_0,x_bar,h,T,M,A,sigma,W,tol)
% X = square_bernoulli(X_0,h,T,M,A,sigma,W). 

N = ceil(T/h + 1);

f = @(x) -A*(x-x_bar);
g = @(x) sigma;
tau = T * ones(M,1);
out = zeros(M,1);

for j = 1:M
    x_old = X_0;
    w = W(2*j-1:2*j,:);
    
    for i = 1:N-1
        x_new = EM_one_step(x_old,f,g,w(:,i+1)-w(:,i),h);
      
        if abs(x_new(1)) > 1 || abs(x_new(2)) > 1
            tau(j) = h*i;
            out(j) = 1;
            break
        else
            [r_1,r_2,r_3,r_4] = compute_distance(x_new);
            bool = critical_distance([r_1,r_2,r_3,r_4],tol);
            p = compute_exit_probability(x_old,x_new,bool,sigma,h);
            if  isempty(find(p > 0.5,1)) == 0
                tau(j) = h*i;
                out(j) = 1;
                break
            end
        end
        x_old = x_new;
    end
end

exp_out = mean(out);
exp_tau = mean(tau);

end