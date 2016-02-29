function p = compute_exit_probability(X,bool,sigma,h)

X_1 = X(:,1);
X_2 = X(:,2);
% count walls 1 left, 2 bottom, 3 right, 4 up

f_p = @(n,z_1,z_2,x_1,x_2) exp(-2 * (n' * (x_1 - z_1)) * (n' * (x_2 - z_1)) / (sigma^2 * h * n' * n));

n_1 = [-1;0];
z_1 = [-1;X_1(2)];
z_2 = [-1;X_2(2)];
p(1) = bool(1) * f_p(n_1,z_1,z_2,X_1,X_2);

n_2 = [0;-1];
z_1 = [X_1(1);-1];
z_2 = [X_2(1);-1];
p(2) = bool(2) * f_p(n_2,z_1,z_2,X_1,X_2);

n_1 = [1;0];
z_1 = [1;X_1(2)];
z_2 = [1;X_2(2)];
p(3) = bool(3) * f_p(n_1,z_1,z_2,X_1,X_2);

n_2 = [0;1];
z_1 = [X_1(1);1];
z_2 = [X_2(1);1];
p(4) = bool(2) * f_p(n_2,z_1,z_2,X_1,X_2);


end