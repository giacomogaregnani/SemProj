A = realizationRF(5,0.05,0.5,1,1);
NGridA = sqrt(length(A));
sizeA = NGridA;
A = reshape(A,NGridA,NGridA);
deltaA = 2 / (NGridA - 1);

dlmwrite('Matrix.txt',sizeA)
dlmwrite('Matrix.txt',deltaA,'-append')
dlmwrite('Matrix.txt',A,'delimiter','\t','precision',3,'-append')