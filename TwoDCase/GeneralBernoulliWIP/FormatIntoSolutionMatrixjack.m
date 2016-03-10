function [U] = FormatIntoSolutionMatrixjack(U_vec,N_X,N_Y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
U = zeros(N_Y,N_X);
for i = 2:N_Y - 2
   
    U(i,2:end-1) = U_vec((i-1)*(N_X-2)+(1:(N_X - 2)) );
  
end


