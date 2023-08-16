clear all;
close all;
clc;

cc = hsv(3);
baseName = 'xAverage_variance_tSteps1000_deltat0.000050_seed';
for seed = 1:3
        fileName = [baseName, num2str(seed), '.txt'];
        data = load(fileName);
        numberOfsteps = data(:,1);
        xAverage = data(:,2);
        variance = data(:,3);  
        figure(1)
        hold on
        plot(numberOfsteps, xAverage, 'color', cc(seed,:))
        xlabel(['Time step t'],'interpreter','latex','FontSize',13)
        ylabel(['$\langle x(t) \rangle$'],'interpreter','latex','FontSize',13)
        figure(2)
        hold on
        plot(numberOfsteps, variance, 'color', cc(seed,:))
        xlabel(['Time step t'],'interpreter','latex','FontSize',13)
        ylabel(['$\sigma^2$'],'interpreter','latex','FontSize',13)
end
figure(1)
%axis([0. 200. 0. 0.35])
legend('1. seed', '2. seed', '3. seed')
figure(2)
%axis([0. 200. 0. 0.08])
legend('1. seed', '2. seed', '3. seed')
