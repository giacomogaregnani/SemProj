function [Ux, Uy, deltaA] = SolveDarcyFF(A,pInlet,plotfields)

system('> SolFiles/*')

sizeA = size(A,1);
deltaA = 2 / (sizeA - 1);

dlmwrite('SolFiles/Matrix.txt',sizeA)
dlmwrite('SolFiles/Matrix.txt',deltaA,'-append')
dlmwrite('SolFiles/Matrix.txt',A,'delimiter','\t','precision',3,'-append')
dlmwrite('SolFiles/Matrix.txt',pInlet,'delimiter','\t','precision',3,'-append')

system('FreeFem++ Darcy.edp');

Ux = dlmread('SolFiles/Ux.txt');
Uy = dlmread('SolFiles/Uy.txt');

X = -1 + deltaA/2 : deltaA : 1 - deltaA/2;
[XX,YY] = meshgrid(X,X);

if plotfields == 'True'
    figure
    surf(XX,YY,Ux,'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
    figure
    surf(XX,YY,Uy,'EdgeColor','none')
    xlabel('x')
    ylabel('y')
end
