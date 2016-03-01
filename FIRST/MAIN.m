% verify the convergence of discrete EM method for killed diffusion
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
h = 15./2.^[0:4];
M = 100;
n_iter = length(h);
exp_out_d = zeros(n_iter,1);
exp_tau_d = exp_out_d;
exp_out_c = exp_out_d;
exp_tau_c = exp_out_d;

% compute the brownian motion
h_ref = h(end)/10;
W = brownian_motion_2D(0,T,h_ref,M);

for j = 1:n_iter
    [exp_out_d(j),exp_tau_d(j)] = square_naive_memory(X_0,x_bar,h(j),T,M,A,sigma,W(:,1:h(j)/h_ref:end));
    [exp_out_c(j),exp_tau_(j)] = square_bernoulli(X_0,x_bar,h(j),T,M,A,sigma,W(:,1:h(j)/h_ref:end),norm(A*h(j)));
end

[X_ex,exp_out_exact,exp_tau_exact] = exact_expectation(A,sigma,X_0,x_bar,T,h_ref,M,W);

err_out_c = abs(exp_out_d - exp_out_exact);
err_out_d = abs(exp_out_c - exp_out_exact);
err_tau_c = abs(exp_tau_d - exp_tau_exact);
err_tau_d = abs(exp_tau_c - exp_tau_exact);


% Plot convergence
figure
loglog(h,err_out_c,'o-')
hold on
loglog(h,sqrt(h),'--')
loglog(h,err_out_d,'<-')
loglog(h,h,'>-')
h_legend=legend('err^{\phi}_d','h^{0.5}','err^{\phi}_c','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('log(h)')
grid on

figure
loglog(h,err_tau_c,'o-')
hold on
loglog(h,sqrt(h),'--')
loglog(h,err_out_d,'<-')
loglog(h,h,'>-')
h_legend=legend('err^{\tau}_d','h^{0.5}','err^{\tau}_c','h');
set(h_legend,'Location','northwest','FontSize',13);
xlabel('log(h)')
grid on

% Compute orders of convergence
orders_out = log2(err_out_c(1:end-1)./err_out_c(2:end));
orders_tau = log2(err_tau_c(1:end-1)./err_tau_c(2:end));
orders_out_b = log2(err_out_d(1:end-1)./err_out_d(2:end));
orders_tau_b = log2(err_tau_d(1:end-1)./err_tau_d(2:end));

