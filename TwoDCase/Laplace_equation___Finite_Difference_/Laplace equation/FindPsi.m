function [phi] = FindPsi(X,Y,phi,A,B,C,D,obstruction_Dimension)


[inside] = CheckInsideObstruction(X,Y,obstruction_Dimension);
Error = 0;
for count = 1:1000
    phi_prev = phi;
    if obstruction_Dimension(1,1) == 4
        for r = 2:(size(X,1)-1)
            for c = 2:(size(X,2)-1)
                value = FindPsi_subfn(phi,r,c,A,B,C,D);phi(r,c) = value;
            end
        end
    elseif obstruction_Dimension(1,1) == 1
        for r = 2:(size(X,1)-1)
            for c = 2:(size(X,2)-1)
                if inside(r,c) == 0
                    value = FindPsi_subfn(phi,r,c,A,B,C,D);phi(r,c) = value;
                end
            end
        end
    end
    error = abs(phi - phi_prev);
    Error(count) = sum(sum(error))/numel(error);
    if Error(count) <= 0.00001;break;end
end



function [value] = FindPsi_subfn(phi,r,c,A,B,C,D)

term1 = phi(r,c-1)/A(r,c)/(A(r,c)+B(r,c));
term2 = phi(r,c+1)/B(r,c)/(A(r,c)+B(r,c));
term3 = phi(r-1,c)/C(r,c)/(C(r,c)+D(r,c));
term4 = phi(r+1,c)/D(r,c)/(C(r,c)+D(r,c));
term5 = (A(r,c)*B(r,c)+C(r,c)*D(r,c))/A(r,c)/B(r,c)/C(r,c)/D(r,c);
value = (term1+term2+term3+term4)/term5;