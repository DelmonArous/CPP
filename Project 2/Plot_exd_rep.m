%omega = 0.01
data = load('data_plot_rep_w001.txt');
rho = data(:,1);
prob_func1 = data(:,2);
prob_func2 = data(:,3);
prob_func3 = data(:,4);

%To get nice plotting
figure()
for i=1:1
    G_state = plot(rho(1:5:200),prob_func1(1:5:200),...
        'ob',rho,prob_func1,'-b');
    hold on
    fst_exci = plot(rho(1:5:200),prob_func2(1:5:200),...
        'or',rho,prob_func2,'-r');
    snd_exci = plot(rho(1:5:200),prob_func3(1:5:200),...
        'og',rho,prob_func3,'-g');

    G_group = hggroup;
    fst_group = hggroup;
    snd_group = hggroup;

    set(G_state,'Parent',G_group)
    set(fst_exci,'Parent',fst_group)
    set(snd_exci,'Parent',snd_group)

    set(get(get(G_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(fst_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(get(get(snd_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Ground state: \lambda_0 = 3','1st excited: \lambda_1 = 7',...
        '2nd excited: \lambda_2 = 11')

    %legend('Ground state: \lambda_0 = 3','','1st excited: \lambda_1 = 7','2nd excited: \lambda_2 = 11')
    xlabel(['$$\rho$$'],'interpreter','latex','FontSize',14)
    ylabel(['$$|\psi|^2$$'],'interpreter','latex','FontSize',14)
end


%omega = 0.5
data = load('data_plot_rep_w050.txt');
rho = data(:,1);
prob_func1 = data(:,2);
prob_func2 = data(:,3);
prob_func3 = data(:,4);

%To get nice plotting
figure()
for i=1:1
    G_state = plot(rho(1:5:200),prob_func1(1:5:200),...
        'ob',rho,prob_func1,'-b');
    hold on
    fst_exci = plot(rho(1:5:200),prob_func2(1:5:200),...
        'or',rho,prob_func2,'-r');
    snd_exci = plot(rho(1:5:200),prob_func3(1:5:200),...
        'og',rho,prob_func3,'-g');

    G_group = hggroup;
    fst_group = hggroup;
    snd_group = hggroup;

    set(G_state,'Parent',G_group)
    set(fst_exci,'Parent',fst_group)
    set(snd_exci,'Parent',snd_group)

    set(get(get(G_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(fst_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(get(get(snd_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Ground state: \lambda_0 = 3','1st excited: \lambda_1 = 7',...
        '2nd excited: \lambda_2 = 11')

    %legend('Ground state: \lambda_0 = 3','','1st excited: \lambda_1 = 7','2nd excited: \lambda_2 = 11')
    xlabel(['$$\rho$$'],'interpreter','latex','FontSize',14)
    ylabel(['$$|\psi|^2$$'],'interpreter','latex','FontSize',14)
end


%omega = 1.0
data = load('data_plot_rep_w100.txt');
rho = data(:,1);
prob_func1 = data(:,2);
prob_func2 = data(:,3);
prob_func3 = data(:,4);

%To get nice plotting
figure()
for i=1:1
    G_state = plot(rho(1:5:200),prob_func1(1:5:200),...
        'ob',rho,prob_func1,'-b');
    hold on
    fst_exci = plot(rho(1:5:200),prob_func2(1:5:200),...
        'or',rho,prob_func2,'-r');
    snd_exci = plot(rho(1:5:200),prob_func3(1:5:200),...
        'og',rho,prob_func3,'-g');

    G_group = hggroup;
    fst_group = hggroup;
    snd_group = hggroup;

    set(G_state,'Parent',G_group)
    set(fst_exci,'Parent',fst_group)
    set(snd_exci,'Parent',snd_group)

    set(get(get(G_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(fst_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(get(get(snd_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Ground state: \lambda_0 = 3','1st excited: \lambda_1 = 7',...
        '2nd excited: \lambda_2 = 11')

    %legend('Ground state: \lambda_0 = 3','','1st excited: \lambda_1 = 7','2nd excited: \lambda_2 = 11')
    xlabel(['$$\rho$$'],'interpreter','latex','FontSize',14)
    ylabel(['$$|\psi|^2$$'],'interpreter','latex','FontSize',14)
end


%omega = 5
data = load('data_plot_rep_w500.txt');
rho = data(:,1);
prob_func1 = data(:,2);
prob_func2 = data(:,3);
prob_func3 = data(:,4);

%To get nice plotting
figure()
for i=1:1
    G_state = plot(rho(1:5:200),prob_func1(1:5:200),...
        'ob',rho,prob_func1,'-b');
    hold on
    fst_exci = plot(rho(1:5:200),prob_func2(1:5:200),...
        'or',rho,prob_func2,'-r');
    snd_exci = plot(rho(1:5:200),prob_func3(1:5:200),...
        'og',rho,prob_func3,'-g');

    G_group = hggroup;
    fst_group = hggroup;
    snd_group = hggroup;

    set(G_state,'Parent',G_group)
    set(fst_exci,'Parent',fst_group)
    set(snd_exci,'Parent',snd_group)

    set(get(get(G_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(fst_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    set(get(get(snd_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Ground state: \lambda_0 = 3','1st excited: \lambda_1 = 7',...
        '2nd excited: \lambda_2 = 11')

    %legend('Ground state: \lambda_0 = 3','','1st excited: \lambda_1 = 7','2nd excited: \lambda_2 = 11')
    xlabel(['$$\rho$$'],'interpreter','latex','FontSize',14)
    ylabel(['$$|\psi|^2$$'],'interpreter','latex','FontSize',14)
end
