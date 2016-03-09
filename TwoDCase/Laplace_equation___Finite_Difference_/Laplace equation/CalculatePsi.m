function [StreamFunction] = CalculatePsi(X,Y,StreamFunction,FM,obstruction_Dimension)

rLimit = size(X,1)-1;cLimit = size(X,2)-1;
for count = 1:1000
    if mod(count,250)==0
        fprintf('Current iteration no. is >> %d \n',count);
    end
    for c = 2:cLimit
        for r = 2:rLimit
            if obstruction_Dimension(1) == 1
                if FM(r,c) == 1
                    a = abs(X(r,c) - X(r,c-1));
                    b = abs(X(r,c+1) - X(r,c));
                    cc= abs(Y(r,c) - Y(r-1,c));
                    d = abs(Y(r+1,c) - Y(r,c));
                    t1 = StreamFunction(r,c-1)/(a*(a+b));
                    t2 = StreamFunction(r,c+1)/(b*(a+b));
                    t3 = StreamFunction(r-1,c)/(cc*(cc+d));
                    t4 = StreamFunction(r+1,c)/(d*(cc+d));
                    t5 = 1/(a*b) + 1/(cc*d);
                    StreamFunction(r,c) = (t1 + t2 + t3 + t4)/t5;
                end
            elseif obstruction_Dimension(1) == 4
                a = abs(X(r,c) - X(r,c-1));
                b = abs(X(r,c+1) - X(r,c));
                cc= abs(Y(r,c) - Y(r-1,c));
                d = abs(Y(r+1,c) - Y(r,c));
                t1 = StreamFunction(r,c-1)/(a*(a+b));
                t2 = StreamFunction(r,c+1)/(b*(a+b));
                t3 = StreamFunction(r-1,c)/(cc*(cc+d));
                t4 = StreamFunction(r+1,c)/(d*(cc+d));
                t5 = 1/(a*b) + 1/(cc*d);
                StreamFunction(r,c) = (t1 + t2 + t3 + t4)/t5;
            end
        end
    end
    if obstruction_Dimension(1) == 1
        cEnd = cLimit + 1;
        for r = 2:rLimit
            a = abs(X(r,cEnd) - X(r,cEnd-1));
            cc= abs(Y(r,cEnd) - X(r-1,cEnd));
            d = abs(Y(r+1,cEnd) - Y(r,cEnd));
            t1 = StreamFunction(r,cEnd-1)/(a*(a+a));
            t3 = StreamFunction(r-1,cEnd)/(cc*(cc+d));
            t4 = StreamFunction(r+1,cEnd)/(d*(cc+d));
            t5 = 1/(a*a) + 1/(cc*d);
            StreamFunction(r,cEnd) = ( 2*t1 + t3 + t4 ) / t5;
        end
    elseif obstruction_Dimension(1) == 4
        a = abs(X(r,c) - X(r,c-1));
        b = abs(X(r,c+1) - X(r,c));
        cc= abs(Y(r,c) - Y(r-1,c));
        d = abs(Y(r+1,c) - Y(r,c));
        t1 = StreamFunction(r,c-1)/(a*(a+b));
        t2 = StreamFunction(r,c+1)/(b*(a+b));
        t3 = StreamFunction(r-1,c)/(cc*(cc+d));
        t4 = StreamFunction(r+1,c)/(d*(cc+d));
        t5 = 1/(a*b) + 1/(cc*d);
        StreamFunction(r,c) = (t1 + t2 + t3 + t4)/t5;
    end
end