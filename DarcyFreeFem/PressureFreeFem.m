CorLen = 0.05;
LMax = 6;
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