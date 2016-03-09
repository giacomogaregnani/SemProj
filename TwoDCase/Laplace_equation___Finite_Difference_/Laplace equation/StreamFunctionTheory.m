x = 0.5:0.5:6;
y = 0.5:0.5:4;
h = 0;
k = 0;
[X,Y] = meshgrid(x,y);
U = 1;
R = 1;
rr = zeros(size(X,1),size(X,2));
theta = zeros(size(X,1),size(X,2));
psi=zeros(size(X,1),size(X,2));
for r = 1:size(X,1)
    for c = 1:size(X,2)
        rr(r,c) = ((h-X(r,c))^2+(k-Y(r,c))^2)^0.5;
%         theta(r,c) = atan((k-Y(r,c))/(h-X(r,c)));
        if X(r,c)>=h && Y(r,c)>=k
            theta(r,c) = 0+atan((k-Y(r,c))/(h-X(r,c)));
        elseif X(r,c)<h && Y(r,c)>k
            theta(r,c) = 90+atan((k-Y(r,c))/(h-X(r,c)));
        elseif X(r,c)<h && Y(r,c)<k
            theta(r,c) = 180+atan((k-Y(r,c))/(h-X(r,c)));
        elseif X(r,c)>h && Y(r,c)<k    
            theta(r,c) = 270+atan((k-Y(r,c))/(h-X(r,c)));
        end 
        psi(r,c) = U*(1-(R/rr(r,c))^2)*rr(r,c)*sind(theta(r,c));
    end
end
figure
contourf(X,Y,psi,25);colorbar;hold on
axis square