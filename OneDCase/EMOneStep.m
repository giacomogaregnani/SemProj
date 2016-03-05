function X = EMOneStep(X,f,g,dW,h)

X = X + f(X)*h + g(X)*dW;

end