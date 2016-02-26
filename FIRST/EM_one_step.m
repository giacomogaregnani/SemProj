function X = EM_one_step(X,f,g,dW,h)

X = X + f(X)*h + g(X)*dW;

end