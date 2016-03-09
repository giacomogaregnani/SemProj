function [ChannelWall,obst_coord] = CreateChannelFLowGeometry(channel_Dimension,obstruction_Dimension)

ChannelWall = [0                      +channel_Dimension(3,1)/2;
               0                      +channel_Dimension(2,1)/2;
               channel_Dimension(1,1) +channel_Dimension(2,1)/2;
               0                      -channel_Dimension(3,1)/2;
               0                      -channel_Dimension(2,1)/2;
               channel_Dimension(1,1) -channel_Dimension(2,1)/2];
obst_coord = 0;
hold on,box on,axis square
plot(ChannelWall(1:3,1),ChannelWall(1:3,2),'r','LineWidth',3)
plot(ChannelWall(1:3,1),ChannelWall(1:3,2),'ko','MarkerSize',8,'MarkerFaceColor',[0 0 0])
plot(ChannelWall(4:6,1),ChannelWall(4:6,2),'r','LineWidth',3)
plot(ChannelWall(4:6,1),ChannelWall(4:6,2),'ko','MarkerSize',8,'MarkerFaceColor',[0 0 0])
if obstruction_Dimension(1,1) == 1
    h = obstruction_Dimension(2,1);k = obstruction_Dimension(6,1);
    X = h + (obstruction_Dimension(3,1)/2)*cosd(0:1:360);
    Y = k + (obstruction_Dimension(3,1)/2)*sind(0:1:360);obst_coord = [X ; Y];
    plot(X,Y,'k','LineWidth',3);fill(X,Y,'g');pause(0.001);axis equal
elseif obstruction_Dimension(1,1) == 2
    X = (obstruction_Dimension(2,1)/2)*cosd(0:1:360);
    Y = (obstruction_Dimension(3,1)/2)*sind(0:1:360);XY = zeros(numel(X),2);
    for c = 1 : numel(X)
       XY(c,:) = [X(c) Y(c)]*[+cosd(obstruction_Dimension(5,1)) sind(obstruction_Dimension(5,1));
                              -sind(obstruction_Dimension(5,1)) cosd(obstruction_Dimension(5,1))];
    end
    X = XY(:,1) + obstruction_Dimension(2,1);Y = XY(:,2) + obstruction_Dimension(6,1);
    obst_coord = [X ; Y];clear XY
    plot(X,Y,'k','LineWidth',3);fill(X,Y,'g');pause(0.001);axis equal
elseif obstruction_Dimension(1,1) == 4 
end