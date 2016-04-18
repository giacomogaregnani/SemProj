CorLen = 0.05;
LMax = 8;
sigmaA = 1;
nu = 0.5;

A = realizationRF(LMax,CorLen,nu,sigmaA,1);
NGridA = sqrt(length(A));
sizeA = NGridA;
A = reshape(A,NGridA,NGridA);
deltaA = 2 / (NGridA - 1);
pIn = 1;

dlmwrite('Matrix.txt',sizeA)
dlmwrite('Matrix.txt',deltaA,'-append')
dlmwrite('Matrix.txt',A,'delimiter','\t','precision',3,'-append')
dlmwrite('Matrix.txt',pIn,'delimiter','\t','precision',3,'-append')

system('FreeFem++ Darcy.edp');

Ux = dlmread('Ux.txt');
Uy = dlmread('Uy.txt');

X = -1 + deltaA/2 : deltaA : 1 - deltaA/2;
[XX,YY] = meshgrid(X,X);

surf(XX,YY,Ux,'EdgeColor','none')
xlabel('x')
ylabel('y')

figure
surf(XX,YY,Uy,'EdgeColor','none')
xlabel('x')
ylabel('y')
