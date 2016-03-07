% histogram parameters
numOfBins = 50;

WEnd = W(:,end);

% normalized histogram plot 

vect = -15:0.01:15;
y = pdf('norm',vect,0,Time(2)^0.5);


[histFreq, histXout] = hist(WEnd, numOfBins);
binWidth = histXout(2)-histXout(1);
bar(histXout, histFreq/binWidth/sum(histFreq));
hold on
plot(vect,y,'k','LineWidth',2)

