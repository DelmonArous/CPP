clear all;
close all;
clc;

data_rk4 = load('data_3body_m1000.txt');
%data_energy = load('energy_esc_n1000.txt');
%data_angmom = load('angmom_esc_n1000.txt');
% t = data_rk4(:,1);
x_sun = data_rk4(:,1);
y_sun = data_rk4(:,2);
x_earth = data_rk4(:,3);
y_earth = data_rk4(:,4);
x_jupiter = data_rk4(:,5);
y_jupiter = data_rk4(:,6);

% x_sun = data_rk4(:,1);
% y_sun = data_rk4(:,2);
% x_mercury = data_rk4(:,3);
% y_mercury = data_rk4(:,4);
% x_venus = data_rk4(:,5);
% y_venus = data_rk4(:,6);
% x_earth = data_rk4(:,7);
% y_earth = data_rk4(:,8);
% x_mars = data_rk4(:,9);
% y_mars = data_rk4(:,10);
% x_jupiter = data_rk4(:,11);
% y_jupiter = data_rk4(:,12);
% x_saturn = data_rk4(:,13);
% y_saturn = data_rk4(:,14);
% x_uranus = data_rk4(:,15);
% y_uranus = data_rk4(:,16);
% x_neptune = data_rk4(:,17);
% y_neptune = data_rk4(:,18);
% x_pluto = data_rk4(:,19);
% y_pluto = data_rk4(:,20);

% E_p = data_energy(:,2);
% E_k = data_energy(:,3);
% E_tot = data_energy(:,4);
% L = data_angmom(:,2);

% To get nice plotting
% figure()
% plot(t, x, '-b');
% hold on
% plot(t, y, '-r');
% h = legend('x','y');
% xlabel(['$\hat{t}\,[yr]$'],'interpreter','latex','FontSize',13)
% ylabel(['Distance from sun [AU]'],'interpreter','latex','FontSize',13)
% 
% % Find the first text object (x) and change it.
% h1 = findobj(get(h,'Children'),'String','x');
% set(h1,'String','$\hat{x}$','Interpreter','latex')
% % Find the second text object (y) and change it.
% h2 = findobj(get(h,'Children'),'String','y');
% set(h2,'String','$\hat{y}$','Interpreter','latex')

figure()
plot(x_sun, y_sun, '-r', x_earth, y_earth, '-b', x_jupiter, y_jupiter, '-k')
legend('Sun', 'Earth', 'Jupiter')
% plot(x_sun, y_sun, '-r', x_mercury, y_mercury, '-b',...
%     x_venus, y_venus, '-g', x_earth, y_earth, '-y',...
%     x_mars, y_mars, '-k', x_jupiter, y_jupiter, '-c',...
%     x_saturn, y_saturn, '-m', x_uranus, y_uranus, '--y',...
%     x_neptune, y_neptune, '--k', x_pluto, y_pluto, '--r')
% legend('Sun', 'Mercury', 'Venus', 'Earth', 'Mars',...
%     'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto')
xlabel(['$x\,[AU]$'],'interpreter','latex','FontSize',13)
ylabel(['$y\,[AU]$'],'interpreter','latex','FontSize',13)

% figure()
% E_0 = E_tot(1);
% E_p0= E_p(1);
% E_k0= E_k(1);
% plot(t, log10(abs((E_tot-E_0)/E_0)), '-b')
% legend('Total energy')
% xlabel(['$\hat{t}\,[yr]$'],'interpreter','latex','FontSize',13)
% ylabel(['$\log_{10}(\epsilon)$'],'interpreter','latex','FontSize',13)
% 
% figure()
% L_0 = L(1);
% plot(t, log10(abs((L-L_0)/L_0)), '-b')
% legend('Angular momentum')
% xlabel(['$\hat{t}\,[yr]$'],'interpreter','latex','FontSize',13)
% ylabel(['$\log_{10}(\epsilon)$'],'interpreter','latex','FontSize',13)


