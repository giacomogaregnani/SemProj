function X = EMOneStepDarcy(X,u,sigma,dW,h)

X = X + u*h + sigma*dW;

end