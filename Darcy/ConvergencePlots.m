function ConvergencePlots(errBernoulliTau, errNaiveTau, tauNaive, tauBernoulli, delta, N, Time)

% Plot Values of Tau
h = Time(2)./N;
semilogx(h, tauNaive(1,:),'o--')
legend('Coarse', 'Medium', 'Fine')
xlabel('h')
ylabel('\tau')
grid on

figure
semilogx(h, tauBernoulli(1,:),'o--')
hold on
semilogx(h, tauBernoulli(2,:),'*--')
semilogx(h, tauBernoulli(3,:),'<--')
legend('Coarse', 'Medium', 'Fine', 'Location', 'NW')
xlabel('h')
ylabel('\tau')
grid on

% Plot Errors
% wrt deltaA

figure
loglog(delta(1:end-1), errBernoulliTau(1:end-1, end), '*--')

% wrt h
figure
loglog(h(1:end-1), errBernoulliTau(end, 1:end-1), 'o--')
hold on
loglog(h(1:end-1), errNaiveTau(end, 1:end-1), 'o--')

end