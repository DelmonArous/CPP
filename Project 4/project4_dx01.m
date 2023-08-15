close all;
clear all;
clc;

%sette inn riktig n for � f� det tidspunktet vi vil ha
for n = [42,202]
    %Analytic plots
    analytic = load('analytic_dx01_dt0001.txt');
    x = analytic(1,:);
    v_analytic_t1 = analytic(n,:);

    %The different methods
    file_name1 = 'explicit_dx01_dt0001.txt';
    scheme1 = load(file_name1);
    file_name2 = 'implicit_dx01_dt0001.txt';
    scheme2 = load(file_name2);
    file_name3 = 'CN_dx01_dt0001.txt';        
    scheme3 = load(file_name3);

    v_scheme1_t1 = scheme1(n,:);
    v_scheme2_t1 = scheme2(n,:);
    v_scheme3_t1 = scheme3(n,:);

    figure()
    plot(x, v_analytic_t1, '-ok')
    xlabel('$x$','interpreter','latex','FontSize',13)
    ylabel('$u(x,t)$','interpreter','latex','FontSize',13)
    axis([0 1 0 1.9]);
    hold on
    plot(x, v_scheme1_t1, '+r')
    plot(x, v_scheme2_t1, '*g')
    plot(x, v_scheme3_t1, 'xb')
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
