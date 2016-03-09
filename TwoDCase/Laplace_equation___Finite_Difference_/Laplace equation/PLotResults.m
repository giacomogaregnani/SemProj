function [] = PLotResults(X,Y,StreamFunction,channel_Dimension,obstruction_Dimension,obst_coord)

figure 
subplot(2,1,1),surf(X,Y,StreamFunction),grid on,box on,axis equal
xlabel('Channel length --- >');ylabel('Channel width --- >')
title('Stream function, \psi'),pause(0.001)

%%
subplot(2,1,2),surf(X,Y,gradient(StreamFunction)),grid on,box on,axis equal
xlabel('Channel length --- >');ylabel('Channel width --- >')
title('Gradient of Stream Function, \nabla\psi'),pause(0.001)

figure
mesh(X,Y,StreamFunction),hold on
plot3(X,Y,StreamFunction,'k.','MarkerSize',12)
grid on,box on,axis equal

if obstruction_Dimension(1,1) == 4
    figure
    subplot(2,1,1),contour(X,Y,StreamFunction,100)
    xlabel('Channel length --- >');ylabel('Channel width --- >')
    title('Stream function, \psi'),pause(0.001),axis equal
    SetAxis(channel_Dimension);
    
    subplot(2,1,2),contour(X,Y,gradient(StreamFunction),100)
    xlabel('Channel length --- >');ylabel('Channel width --- >')
    title('Gradient of Stream Function, \nabla\psi'),pause(0.001),axis equal
    SetAxis(channel_Dimension);
    
elseif obstruction_Dimension(1,1) == 1
    figure
    subplot(2,1,1),contour(X,Y,StreamFunction,2000),hold on
    plot(obst_coord(1,:),obst_coord(2,:),'k-','LineWidth',3),fill(X,Y,'g'),alpha(0.5);
    xlabel('Channel length --- >');ylabel('Channel width --- >')
    title('Stream function, \psi'),pause(0.001),axis equal
    SetAxis(channel_Dimension);
    
    subplot(2,1,2),contour(X,Y,gradient(StreamFunction),2000),hold on
    plot(obst_coord(1,:),obst_coord(2,:),'k-','LineWidth',3),fill(X,Y,'g'),alpha(0.5);
    xlabel('Channel length --- >');ylabel('Channel width --- >')
    title('Gradient of Stream Function, \nabla\psi'),pause(0.001),axis equal
    SetAxis(channel_Dimension);
end

figure
contour(X,Y,StreamFunction,2000)

figure
contour(X,Y,gradient(StreamFunction),2000)

function SetAxis(channel_Dimension)
xmin = 0;
xmax = channel_Dimension(1,1);
ymin = -channel_Dimension(2,1)/2;
ymax = channel_Dimension(2,1)/2;
box on,axis equal
axis([xmin xmax ymin ymax])
