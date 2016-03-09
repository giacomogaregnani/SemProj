contour(X,Y,StreamFunction,800)
axis equal
xlabel('Channel Length')
ylabel('Channel Width')
title('Value of Stream Function \psi')
colorbar
print('-djpeg100','6a_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_contour.jpeg')
%%
figure
start = floor(size(X,2)/1.75);
NoOfContours = 100;
contour(X(:,start:size(X,2)),Y(:,start:size(X,2)),StreamFunction(:,start:size(X,2)),NoOfContours)
axis equal
xlabel('Channel Length')
ylabel('Channel Width')
title('Value of Stream Function \psi')
colorbar
print('-djpeg100','6b_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_contour.jpeg')
%%
figure
startx = floor(size(X,1)/1.25);
starty = 1;
endy = floor(size(X,2)/6);
NoOfContours = 100;
contour(X(startx:size(X,1),1:endy),Y(startx:size(X,1),1:endy),StreamFunction(startx:size(X,1),1:endy),NoOfContours)
axis equal
xlabel('Channel Length')
ylabel('Channel Width')
title('Stream Function \psi around the corner')
colorbar
print('-djpeg100','6f_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_contour.jpeg')