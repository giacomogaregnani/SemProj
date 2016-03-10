clc
clear

close all

X_min = -1;
X_max = 1;
N_X = 100;
X = linspace(X_min,X_max,N_X);
Y_min = -1;
Y_max = 1;
N_Y = 100;
Y = linspace(Y_min,Y_max,N_Y);
g = [1 0; 0 1];

[A,F] = AssembleSpatialOperatorWithBC(N_X,N_Y,X,Y,g);
u_vec = A\F;
%formats into matrix and adds BC
%[U] = FormatIntoSolutionMatrixjack(u_vec,N_X,N_Y);
U = zeros(N_X,N_Y);
U(2:end-1,2:end-1) = reshape(u_vec,N_X-2,N_Y-2);
[XX,YY] = meshgrid(X,Y);
surf(XX,YY,U)

a = (N_X - 1)^2 / 4 * gallery('poisson',N_X-2);
u = -a\F;
U2 = zeros(N_X,N_Y);
U2(2:end-1,2:end-1) = reshape(u,N_X-2,N_Y-2);

figure
surf(XX,YY,U - U2)

