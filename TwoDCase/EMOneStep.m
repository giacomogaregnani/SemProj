function X = EMOneStep(X,f,g,dW,h)

X = X + f(X(1),X(2))*h + g(X(1),X(2))*dW;

end