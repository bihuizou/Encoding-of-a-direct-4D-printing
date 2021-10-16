clc
clear all
close all
load('ErrorData_2_3_T85_theta0_theta45.mat')
a = 0.01:0.001:0.04;
b = 0.01:0.001:0.04;

figure(1)
frontsize = 15;
subplot(1,3,1)
pcolor(a,b,theta_error_1');
colorbar
shading flat
set(gca,'YaxisLocation','right','Box','off','FontSize',frontsize,'FontWeight','bold','FontName','Calibri');
xlabel('Length a [m]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')
ylabel('Length b [m]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')



subplot(1,3,2)
pcolor(a,b,a_error_1');
colorbar
shading flat
set(gca,'YaxisLocation','right','Box','off','FontSize',frontsize,'FontWeight','bold','FontName','Calibri');
xlabel('Length a [m]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')
ylabel('Length b [m]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')

subplot(1,3,3)
pcolor(a,b,b_error_1');
colorbar
shading flat
set(gca,'YaxisLocation','right','Box','off','FontSize',frontsize,'FontWeight','bold','FontName','Calibri');
xlabel('Length a [m]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')
ylabel('Length b [m]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')