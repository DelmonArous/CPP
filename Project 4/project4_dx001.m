close all;
clear all;
clc;

%sette inn riktig n for å få det tidspunktet vi vil ha
for n = [4002,20002]
    %Analytic plots
    analytic = load('analytic_dx001_dt000001.txt');
    x = analytic(1,:);
    v_analytic_t1 = analytic(n,:);

    %The different methods
    file_name1 = 'explicit_dx001_dt000001.txt';
    scheme1 = load(file_name1);
    v_scheme1_t1 = scheme1(n,:);
    file_name1 = 'implicit_dx001_dt000001.txt';
    scheme2 = load(file_name1);
    v_scheme2_t1 = scheme2(n,:);
    file_name1 = 'CN_dx001_dt000001.txt';        
    scheme3 = load(file_name1);
    v_scheme3_t1 = scheme3(n,:);
    

    figure()
    plot(x(1:5:length(x)), v_analytic_t1(1:5:length(x)), '-ok')
    xlabel('$x$','interpreter','latex','FontSize',13)
    ylabel('$u(x,t)$','interpreter','latex','FontSize',13)
    axis([0 1 0 1.9]);
    hold on
    plot(x(1:5:length(x)), v_scheme1_t1(1:5:length(x)), '+r')
    plot(x(1:5:length(x)), v_scheme2_t1(1:5:length(x)), '*g')
    plot(x(1:5:length(x)), v_scheme3_t1(1:5:length(x)), 'xb')
    legend('Analytic solution','Explicit','Implicit','Crank-Nicolson');

%Regner ut feilen i forhold til analytisk
    for i =1:length(analytic(n,:))
        v_scheme1_t1(i) = abs(v_scheme1_t1(i)-v_analytic_t1(i));
        v_scheme2_t1(i) = abs(v_scheme2_t1(i)-v_analytic_t1(i));
        v_scheme3_t1(i) = abs(v_scheme3_t1(i)-v_analytic_t1(i));
    end
    figure()
    plot(x, v_scheme1_t1, '-r')
    hold on
    plot(x, v_scheme2_t1, '-g')
    plot(x, v_scheme3_t1, '-b')
    legend('Explicit','Implicit','Crank-Nicolson');
    xlabel('$x$','interpreter','latex','FontSize',13)
    ylabel('Absolute error','interpreter','latex','FontSize',13)
    title('Compared to analytic soluton')
end
