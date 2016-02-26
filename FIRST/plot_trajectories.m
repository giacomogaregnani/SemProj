function plot_trajectories(X,M)

plot([-1,-1],[-1,1],'k',[1,1],[-1,1],'k',[-1,1],[-1,-1],'k',[-1,1],[1,1],'k','LineWidth',2)
axis([-1.5 1.5 -1.5 1.5])
hold on

for i = 1:M
    plot(X(2*i-1,:),X(2*i,:))
end

end