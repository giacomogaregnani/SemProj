function [Ux, Uy] = SolveDarcyFF(A,pInlet,plotfields,deltaU)

system('> SolFiles/Ux.txt');
system('> SolFiles/Uy.txt');
system('> SolFiles/Data.txt');

sizeA = size(A,1);
deltaA = 2 / (sizeA - 1);
sizeOut = 2 / deltaU;

dlmwrite('SolFiles/Data.txt', sizeA)
dlmwrite('SolFiles/Data.txt', deltaA, '-append')
dlmwrite('SolFiles/Data.txt', deltaU, '-append')
dlmwrite('SolFiles/Data.txt', A, 'delimiter', '\t', 'precision', 3, '-append')
dlmwrite('SolFiles/Data.txt',pInlet, 'delimiter', '\t', 'precision', 3, '-append')
dlmwrite('SolFiles/Data.txt', sizeOut, '-append')

system('FreeFem++ Darcy.edp');

Ux = dlmread('SolFiles/Ux.txt');
Uy = dlmread('SolFiles/Uy.txt');

X = -1 + deltaU/2 : deltaU : 1 - deltaU/2;
[XX,YY] = meshgrid(X,X);

if strcmp(plotfields, 'True') == 1
    figure
    surf(XX,YY,Ux,'EdgeColor','none')
    xlabel('x')
    ylabel('y')
    
    figure
    surf(XX,YY,Uy,'EdgeColor','none')
    xlabel('x')
    ylabel('y')
end
