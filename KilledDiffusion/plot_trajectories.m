function plot_trajectories(M,h,X_0,x_bar,A,sigma)

figure
p(1) = plot([-1,-1],[-1,1],'k','LineWidth',2)
hold on
plot([1,1],[-1,1],'k',[-1,1],[-1,-1],'k',[-1,1],[1,1],'k','LineWidth',2)

axis([-1.2 1.2 -1.2 1.2])
axis equal
hold on

W = brownian_motion_2D(0,30,h(end),M);

for i = 1:length(h)
    X = exact_expectation(A,sigma,X_0,x_bar,30,h(i),M,W(:,1:h(i)/h(end):end));
    for j = 1:M
        if i == 1
            p(i+1) = plot(X(2*j-1,:),X(2*j,:),'LineWidth',2)
        else
            p(i+1) = plot(X(2*j-1,:),X(2*j,:))
    end
    
end

legend(p,{'Domain','big timestep','small timestep'})
end