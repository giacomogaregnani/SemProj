% verify the convergence of corrected EM method for killed diffusion
% problems. 

clear
close all
clc

% Solution of the transport diffusion SDE with velocity field v = -Ax;

% define the problem
A = 0.01 * [1,0;0,3];
sigma = 0.1;
X_0 = 0.5 * ones(2,1);
x_bar = [-0.9;-0.9];

% define the integration
T = 30;
h = 15./2.^[0:5];
M = 2000;
n_iter = length(h);
exp_out = zeros(n_iter,1);
exp_tau = exp_out;

% compute the brownian motion
h_ref = h(end)/10;
W = brownian_motion_2D(0,T,h_ref,M);

for j = 1:n_iter
    
    [exp_out(j),exp_tau(j)] = square_bernoulli(X_0,x_bar,h(j),T,M,A,sigma,W(:,1:h(j)/h_ref:end),0.1);
        
end


[X_ex,exp_out_exact,exp_tau_exact] = exact_expectation(A,sigma,X_0,x_bar,T,h_ref,M,W);

err_out = abs(exp_out - exp_out_exact);
err_tau = abs(exp_tau - exp_tau_exact);

figure
loglog(h,err_out,'o-')
hold on
loglog(h,h,'*-')
legend('err_{out}','h','Location','NW')
xlabel('log(h)')
grid on

figure
loglog(h,err_tau,'o-')
hold on
loglog(h,h,'*-')
legend('err_{tau}','h','Location','NW')
xlabel('log(h)')
grid on

orders_out = log2(err_out(1:end-1)./err_out(2:end));
orders_tau = log2(err_tau(1:end-1)./err_tau(2:end));