function InletPoints = IdentifyInletPoints(X,Y,channel_Dimension)

track = 1;
for count = 1:size(Y,1)
    Dist00 = sqrt(X(count,1)^2 + Y(count,1)^2);
    if Dist00 < channel_Dimension(3,1)/2;
        InletPoints(track) = count;
        track = track + 1;
    end
end