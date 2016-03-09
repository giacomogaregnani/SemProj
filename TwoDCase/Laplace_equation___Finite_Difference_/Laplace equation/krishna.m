InitialSettings
%%
[channel_Dimension,obstruction_Dimension,grid_values,InletVel] = INPUT_ChannelFlow();
% tic
%%
[ChannelWall,obst_coord] = CreateChannelFLowGeometry(channel_Dimension,obstruction_Dimension);
%%
[X,Y] = MakeGrid(channel_Dimension,grid_values);
%%
[StreamFunction] = StreamFunction_Define(X,Y,channel_Dimension,InletVel);
%%
[FM] = AllocateFlagValue(X,Y,obstruction_Dimension);
%%
[StreamFunction] = CalculatePsi(X,Y,StreamFunction,FM,obstruction_Dimension);
% PLotResults(X,Y,StreamFunction,channel_Dimension,obstruction_Dimension,obst_coord)
% toc
%%
[FM] = AllocateFlagValue1(X,Y,obstruction_Dimension,grid_values);
%%
[u,v,U] = CalculateVelocity(X,Y,StreamFunction,FM);
% clear obst_coord grid_values h k FM x y

