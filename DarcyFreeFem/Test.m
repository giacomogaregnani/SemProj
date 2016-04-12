A = repmat([1:5],5,1)';
dA = 2 / (size(A,1)-1);

x = 0.9;
y = 1;

xLeft = -1 + (ceil((x+1)/dA)-1)*dA;
xRight = -1 + (ceil(1 + (x+1)/dA)-1)*dA;
yLeft = -1 + (ceil((y+1)/dA)-1)*dA;
yRight = -1 + (ceil(1+(y+1)/dA)-1)*dA;

f = (1/dA^2) * (A(ceil((x+1)/dA),ceil((y+1)/dA))*(xRight-x)*(yRight-y) ...
    + A(ceil(1+(x+1)/dA),ceil((y+1)/dA))*(x-xLeft)*(yRight-y) ...
    + A(ceil((x+1)/dA),ceil(1+(y+1)/dA))*(xRight-x)*(y-yLeft) ...
    + A(ceil(1+(x+1)/dA),ceil(1+(y+1)/dA))*(x-xLeft)*(y-yLeft))
