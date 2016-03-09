function [u,v,U] = CalculateVelocity(X,Y,StreamFunction,FM)

% u = zeros(size(X,1),floor(size(X,2)/4));
% v = zeros(size(X,1),floor(size(X,2)/4));

u = zeros(size(X,1),size(X,2));
v = zeros(size(X,1),size(X,2));
% U = zeros(size(X,1),size(X,2));

%%
for r = 3:size(X,1)-2
    for c = 3:size(X,2)-2
        if FM(r,c) == 1
            u(r,c) =  abs(StreamFunction(r+1,c)-StreamFunction(r-1,c))./abs(Y(r+1,c)-Y(r-1,c));
            v(r,c) = -(StreamFunction(r,c+1)-StreamFunction(r,c-1))./abs(X(r,c+1)-X(r,c-1));
        end
    end
end


% rLimit = size(X,1);
% cLimit = size(X,2);
% 
% r = 3:rLimit-2;c = 3:cLimit-2;
% u(r,c) =  (StreamFunction(r+1,c)-StreamFunction(r-1,c))./(Y(r+1,c)-Y(r-1,c));
% v(r,c) = -(StreamFunction(r,c+1)-StreamFunction(r,c-1))./(X(r,c+1)-X(r,c-1));
%%
U = sqrt(u.^2 + v.^2);
% angle = atan(u./v);
% 
% contour(X(3:(size(X,1)-2),3:(size(X,2)-2)),Y(3:(size(X,1)-2),3:(size(X,2)-2)),u(3:(size(X,1)-2),3:(size(X,2)-2)),200)
% contour(X(3:(size(X,1)-2),3:(size(X,2)-2)),Y(3:(size(X,1)-2),3:(size(X,2)-2)),v(3:(size(X,1)-2),3:(size(X,2)-2)),100)
% contour(X(3:(size(X,1)-2),3:(size(X,2)-2)),Y(3:(size(X,1)-2),3:(size(X,2)-2)),U(3:(size(X,1)-2),3:(size(X,2)-2)),100)