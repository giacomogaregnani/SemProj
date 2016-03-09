xx = X(3:((size(X,1)+1)/2),3:(size(X,2)-2));
yy = Y(3:((size(Y,1)+1)/2),3:(size(Y,2)-2));
UU = U(3:((size(X,1)+1)/2),3:(size(X,2)-2));
uu = u(3:((size(X,1)+1)/2),3:(size(X,2)-2));
vv = v(3:((size(X,1)+1)/2),3:(size(X,2)-2));
%%
for count = 1:size(xx,1)
    plot(xx(count,1:size(xx,2)),UU(count,1:size(xx,2)),xx(count,1:size(xx,2)),UU(count,1:size(xx,2)),'k.')
    hold on
end
axis square
axis tight
grid on
xlabel('Channel Length')
ylabel('Velocity, U')
title('Variation of velocity U along the channel length across the circular cylinder')
print('-djpeg100','6c_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_contour.jpeg')
%%
figure
for count = 1:size(xx,1)
    plot(xx(count,1:size(xx,2)),uu(count,1:size(xx,2)),xx(count,1:size(xx,2)),uu(count,1:size(xx,2)),'k.')
    hold on
end
axis square
axis tight
grid on
xlabel('Channel Length')
ylabel('Velocity, u')
title('Variation of velocity u along the channel length across the circular cylinder')
print('-djpeg100','6d_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_contour.jpeg')
%%
figure
for count = 1:size(xx,1)
    plot(xx(count,1:size(xx,2)),vv(count,1:size(xx,2)),xx(count,1:size(xx,2)),vv(count,1:size(xx,2)),'k.')
    hold on
end
axis square
axis tight
grid on
xlabel('Channel Length')
ylabel('Velocity, v')
title('Variation of velocity v along the channel length across the circular cylinder')
print('-djpeg100','6e_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_contour.jpeg')
%%
figure
for count = 1:size(yy,2)
    plot(yy(1:size(yy,1),count),UU(1:size(yy,1),count),yy(1:size(yy,1),count),UU(1:size(yy,1),count),'k.')
    hold on
end
%%
axis square
axis tight
grid on
xlabel('Channel Width')
ylabel('Velocity, U')
title('Variation of velocity U along the channel width across the circular cylinder')
print('-djpeg100','6g_L6W4I1D3d2O0X0.10Y0.10V1TBpsi1_Lineplot.jpeg')