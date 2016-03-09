function [inside] = CheckInsideObstruction(X,Y,obstruction_Dimension)

inside = zeros(size(X,1),size(X,2));
if obstruction_Dimension == 1
    dist_hk = ((obstruction_Dimension(2,1)-X(:,:)).^2 +...
               (obstruction_Dimension(6,1)-Y(:,:)).^2).^0.5;
    inside(:,:) = dist_hk(:,:)<= obstruction_Dimension(3,1)/2;
elseif obstruction_Dimension == 4
    inside = NaN;
end