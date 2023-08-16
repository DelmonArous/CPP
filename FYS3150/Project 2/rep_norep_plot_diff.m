%omega = 0.01
data = load('data_plot_norep_w500.txt');
rho = data(:,1);
prob_func1_no = data(:,2);
prob_func2_no = data(:,3);
prob_func3_no = data(:,4);

%omega = 0.01
data = load('data_plot_rep_w500.txt');
prob_func1_rep = data(:,2);
prob_func2_rep = data(:,3);
prob_func3_rep = data(:,4);

figure()
plot(rho,prob_func1_rep,'-r',rho,prob_func2_rep,':r')
hold on
plot(rho,prob_func1_no,'-b',rho,prob_func2_no,':b')
legend('Repulsion Ground state','Repulsion 1st excited',...
    'No Repulsion Ground state','No Repulsion 1st excited')
xlabel(['$$\rho$$'],'interpreter','latex','FontSize',14)
ylabel(['$$|\psi|^2$$'],'interpreter','latex','FontSize',14)