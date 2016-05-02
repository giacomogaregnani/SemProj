function ConvergencePlots(errBernoulliTau, errNaiveTau, tauBernoulli, delta, N, Time)

% Plot Values of Tau
h = Time(2)./N;

figure
mrk = 'o*<>sd';
for j = 1 : size(tauBernoulli, 1)
    semilogx(h, tauBernoulli(j,:),'Marker',mrk(j))
    hold on
    legendtext{j} = ['\Delta_u = ', num2str(delta(j))];
end
legend(legendtext, 'Location', 'NW');

% Plot Errors
% wrt Deltau

figure
loglog(delta(1:end-1), errBernoulliTau(1:end-1, end), '*--')
xlabel('\Delta_u')
ylabel('error')

% wrt h
figure
loglog(h(1:end-1), errBernoulliTau(end, 1:end-1), 'o--')
hold on
loglog(h(1:end-1), errNaiveTau(end, 1:end-1), 'o--')

end