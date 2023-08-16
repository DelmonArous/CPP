clear all;
close all;
clc;

data = load('xAverage_variance_gauss_tSteps200_deltat0.005000.txt');
tStep = data(:,1);
xAverage = data(:,2);
variance = data(:,3);
figure()
plot(tStep, xAverage)
xlabel(['Time step t'],'interpreter','latex','FontSize',13)
ylabel(['$\langle x(t) \rangle$'],'interpreter','latex','FontSize',13)
figure()
plot(tStep, variance)
xlabel(['Time step $t$'],'interpreter','latex','FontSize',13)
ylabel(['$\sigma^2$'],'interpreter','latex','FontSize',13)

baseName = 'histogram_gauss_tSteps';
for t = [200]
    fileName = [baseName, num2str(t), '_deltat0.005000.txt'];
    data = load(fileName);
    u = data(:,2);
    x = linspace(0, 1, length(u));
    figure();
    bar(x, u)
    
    legend([num2str(t), ' time steps']);
    xlabel(['Steps $x$'],'interpreter','latex','FontSize',13)
    ylabel(['$u(x,t)$'],'interpreter','latex','FontSize',13)
end