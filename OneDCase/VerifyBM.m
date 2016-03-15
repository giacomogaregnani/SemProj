% Program to verify the correctness of the function BrownianMotion

WEnd = W(:,end);

% normalized histogram plot 

vect = -15:0.01:15;
y = pdf('norm',vect,0,Time(2)^0.5);

numOfBins = 50;
[histFreq, histXout] = hist(WHATYOUWANTTOPLOT, numOfBins);
binWidth = histXout(2)-histXout(1);
bar(histXout, histFreq/binWidth/sum(histFreq));
hold on
plot(vect,y,'k','LineWidth',2)

