function [FM] = AllocateFlagValue(X,Y,obstruction_Dimension)

h = obstruction_Dimension(2);
k = obstruction_Dimension(6);
D =obstruction_Dimension(3);
FM = ones(size(X,1),size(X,2));
FM(:,1) = NaN;
FM(:,size(X,2)) = NaN;
FM(1,:) = NaN;
FM(size(X,1),:) = NaN;

rLimit = size(X,1)-1;
cLimit = size(X,2)-1;
for r = 2:rLimit
    for c = 2:cLimit
        T1 = X(r,c)-h;
        T2 = Y(r,c)-k;
        if ((T1*T1+T2*T2)^0.5) <= (D/2)
            FM(r,c) = 0;
        end
    end
end