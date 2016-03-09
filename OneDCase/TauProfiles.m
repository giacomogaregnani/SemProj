function TauProfiles(V,dV,g,Bounds,BoundCond,W,Time,SpaceSampling)
% Plot the profiles of tau for the exact solution, a DEM and a CEM
% approximation of the solution

InitGuess = linspace(Bounds(1),Bounds(2),SpaceSampling);
f = @(x) -dV(x);
tauNaive = zeros(1,SpaceSampling);
tauBernoulli = tauNaive;
tauExact = tauNaive;

for i = 1 : SpaceSampling
    tauNaive(i) = ComputeExitTimeNaive(InitGuess(i),f,g,Bounds,BoundCond,W,Time);
    tauBernoulli(i) = ComputeExitTimeBernoulli(InitGuess(i),f,g,Bounds,BoundCond,W,Time);
    tauExact(i) = ComputeExitTimeExact(InitGuess(i),V,g,Bounds,BoundCond);
end

figure
hold on
plot(InitGuess,tauNaive,'ro--')
plot(InitGuess,tauBernoulli,'b*--')
plot(InitGuess,tauExact,'k<--')
legend('DEM','CEM','Exact','Location','NW')
grid on
xlabel('X_0')
ylabel('\tau')

end

