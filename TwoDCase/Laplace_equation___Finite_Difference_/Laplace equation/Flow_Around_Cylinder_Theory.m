x = -3:0.2:3;
y = -2:0.2:2;
h = 0;
k = 0;
[X,Y] = meshgrid(x,y);
U = 0.5;
R = 1;

for r = 1:size(X,1)
    for c = 1:size(X,2)
        psi(r,c) = U * (1-(R/rr(r,c))^2)*rr(r,c)*sind(theta(r,c));