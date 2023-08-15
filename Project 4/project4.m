close all;
clear all;
clc;

analytic = load('analytic_dx01_dt0001.txt');
x = analytic(1,:);
v_analytic_t1 = analytic(50,:);
v_analytic_t2 = analytic(5,:);

for i = 1:3
    if (i==1)
        scheme = load('explicit_dx01_dt0001.txt');
    elseif (i==2)
        scheme = load('implicit_dx01_dt0001.txt');
    else
        scheme = load('CN_dx01_dt0001.txt');
    end
    
    v_scheme_t1 = scheme(50,:);
    v_scheme_t2 = scheme(5,:);

    figure();
    plot(x, v_analytic_t1, '-o', x, v_scheme_t1, '*')
    legend('Analytic solution', 'Method');
    xlabel('x','interpreter','latex','FontSize',13)
    ylabel('v','interpreter','latex','FontSize',13)
end