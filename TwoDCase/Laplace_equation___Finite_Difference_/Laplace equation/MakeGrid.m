function [X,Y] = MakeGrid(channel_Dimension,grid_values)

x = 0:grid_values(1):channel_Dimension(1,1);
y = (0:grid_values(2):channel_Dimension(2,1)) - channel_Dimension(2,1)/2;
[X,Y] = meshgrid(x,y);