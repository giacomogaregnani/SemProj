function X = EMOneStep(X,f,sigma,h)

X = X + f(X(1),X(2))*h + sigma*sqrt(h)*randn(2,1);

end