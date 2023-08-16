clear all;
close all;
clc;

analytic = importdata('analytic_t1000_deltat0.000050.txt');
x_analytic = analytic(:,1);
u_analytic = analytic(:,2);

cc = hsv(3);
linS = {'-','--',':'};
baseName = 'histogram_tSteps1000_deltat0.000050_seed';
for seed = 1:3
    fileName = [baseName, num2str(seed), '.txt'];
    data = load(fileName);
    uHist = data(:,2);
    x_data = linspace(0, 1, length(uHist));
    figure(1)
    hold on
    stairs(x_data, uHist, 'color', cc(seed,:))%, 'linestyle', linS{seed})
    xlabel(['Steps $x$'],'interpreter','latex','FontSize',13)
    ylabel(['$u(x,t)$'],'interpreter','latex','FontSize',13)
    figure(2)
    hold on
    plot (x_analytic, abs(u_analytic - uHist), 'color', cc(seed,:))
    xlabel(['Steps $x$'],'interpreter','latex','FontSize',13)
    ylabel(['Absolute error'],'interpreter','latex','FontSize',13)
end
figure(1)
hold on
plot(x_analytic, u_analytic, '-kx')
legend('1. seed', '2. seed', '3. seed', 'analytic');
axis([0. 1. 0. 1.])
title('1000 time steps')
figure(2)
hold on
legend('1. seed', '2. seed', '3. seed');
title('1000 time steps')
