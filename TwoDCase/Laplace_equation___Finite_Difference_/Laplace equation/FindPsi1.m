function [phi,place,insidepoint] = FindPsi1(X,Y,phi,A,B,C,D,obstruction_Dimension)

SzInsider = size(X,1);SzInsidec = size(X,2);
insidepoint = zeros(SzInsider,SzInsidec);

if obstruction_Dimension(1,1) == 1
    place = zeros(size(X,1),size(X,2));
    h = obstruction_Dimension(2,1);k = obstruction_Dimension(6,1);
    for e = 1:numel(X)
        dist = ((X(e)-h)^2 + (Y(e)-k)^2)^0.5;
        if dist <= obstruction_Dimension(3)
            place(e) = e;
        end
    end
    [zr,zc] = find(place == 0,numel(place));
    place(1:(zr(1)-1),:)=[];
    place((zr(numel(zr))+1):SzInsider,:)=[];
    place(:,1:(zc(1)-1))=[];
    place(:,(zc(numel(zc))+1):SzInsidec)=[];
else
    place = 0;
end

Error = 0;
for count = 1:1000
    phi_prev = phi;
    if obstruction_Dimension(1,1) == 1
        for r = 2:(SzInsider-1)
            for c = 2:SzInsidec
                if c ~= SzInsidec
                    ChannelEnd = 0;
                else
                    ChannelEnd = 1;
                end
                
                if place(r,c) == 0 
                    value = FindPsi_subfn(phi,X,r,c,A,B,C,D,ChannelEnd);
                    phi(r,c) = value;
                end
            end
        end
    elseif obstruction_Dimension(1,1) == 4
            for r = 2:(SzInsider-1)
            for c = 2:SzInsidec
                if c ~= SzInsidec;ChannelEnd = 0;else ChannelEnd = 1;end
                value = FindPsi_subfn(phi,X,r,c,A,B,C,D,ChannelEnd);
                phi(r,c) = value;
            end
            end
    end
    error = abs(phi - phi_prev);
    Error(count) = sum(sum(error))/numel(error);
    if Error(count) <= 0.00001;break;end
end

function [value] = FindPsi_subfn(phi,X,r,c,A,B,C,D,ChannelEnd)

if ChannelEnd == 0
    term1 = phi(r,c-1)/A(r,c)/(A(r,c)+B(r,c));
    term2 = phi(r,c+1)/B(r,c)/(A(r,c)+B(r,c));
    term3 = phi(r-1,c)/C(r,c)/(C(r,c)+D(r,c));
    term4 = phi(r+1,c)/D(r,c)/(C(r,c)+D(r,c));
    term5 = ( A(r,c)*B(r,c) + C(r,c)*D(r,c) ) / ( A(r,c)*B(r,c)*C(r,c)*D(r,c) );
    value = (term1+term2+term3+term4)/term5;
elseif ChannelEnd ==1
    term1 = phi(r,c-1)/A(r,c)/(A(r,c)+B(r,c));
    term3 = phi(r-1,c)/C(r,c)/(C(r,c)+D(r,c));
    term4 = phi(r+1,c)/D(r,c)/(C(r,c)+D(r,c));
    term5 = ( A(r,c)*B(r,c) + C(r,c)*D(r,c) ) / ( A(r,c)*B(r,c)*C(r,c)*D(r,c) );
    value = (2*term1+term3+term4)/term5;
end