x = -0.1:0.01:1;y = -1:0.01:0.1;
[X,Y] = meshgrid(x,y);psi = 2*X.*Y;
figure;contourf(X,Y,psi,25);colorbar;hold on
plot([min(x) max(x)],[0 0],'k','LineWidth',3)
plot([0 0],[min(y) max(y)],'k','LineWidth',3)
xlabel('X-axis');ylabel('Y-axis')
title('Stream Function \psi near the corners of walls')
axis square;axis tight;print('-djpeg100','CornerFlow_Theory')