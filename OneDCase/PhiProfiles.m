function PhiProfiles(f,g,Bounds,BoundCond,W,Time,SpaceSampling)
% Plot the profiles of phi for the exact solution, a DEM and a CEM
% approximation of the solution

InitGuess = linspace(Bounds(1),Bounds(2),SpaceSampling);
phiNaive = zeros(1,SpaceSampling);
phiBernoulli = phiNaive;
phiExact = phiNaive;

for i = 1 : SpaceSampling
    [b,phiNaive(i)] = ComputeExitTimeNaive(InitGuess(i),f,g,Bounds,BoundCond,W,Time);
    [a,phiBernoulli(i)] = ComputeExitTimeBernoulli(InitGuess(i),f,g,Bounds,BoundCond,W,Time);
    phiExact(i) = ComputeExitProbExact(InitGuess(i),Time,Bounds,BoundCond,f,g(1));
end

figure
hold on
plot(InitGuess,phiNaive,'ro--')
plot(InitGuess,phiBernoulli,'b*--')
plot(InitGuess,phiExact,'k<--')
legend('DEM','CEM','Exact','Location','NW')
grid on
xlabel('X_0')
ylabel('\tau')

end

