function [StreamFunction] = StreamFunction_Define(X,Y,channel_Dimension,InletVel)

StreamFunction = zeros(size(X,1),size(X,2));
InletPoints = IdentifyInletPoints(X,Y,channel_Dimension);
StreamFunction(InletPoints(1,:)) = InletVel;
% StreamFunction(size(StreamFunction,1),:) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: USE THIS TO DRFINE SURROUNDING BOUNDARY CONDITION
% StreamFunction(1,:) = 1;