function p = compute_exit_probability(X_1,X_2,bool,sigma,h)

% count walls 1 left, 2 bottom, 3 right, 4 up

f_p = @(n,z,x_1,x_2) exp(-2 * (n' * (x_1 - z)) * (n' * (x_2 - z)) / (sigma^2 * h * n' * n));

n = [-1;0];
z = [-1;X_1(2)];
p(1) = bool(1) * f_p(n,z,X_1,X_2);

n = [0;-1];
z = [X_1(1);-1];
p(2) = bool(2) * f_p(n,z,X_1,X_2);

n = [1;0];
z = [1;X_1(2)];
p(3) = bool(3) * f_p(n,z,X_1,X_2);

n = [0;1];
z = [X_1(1);1];
p(4) = bool(2) * f_p(n,z,X_1,X_2);


end