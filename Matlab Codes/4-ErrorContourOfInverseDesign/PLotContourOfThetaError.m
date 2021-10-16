clc
clear all
close all
load('ErrorData_1_2_T85_a30_b10.mat')
t1 = -90:1:90;
t2 = -90:1:90;

figure(1)
frontsize = 18;
subplot(1,3,1)
pcolor(t1,t2,theta_error_1');
colorbar
shading flat
set(gca,'YaxisLocation','right','Box','off','FontSize',frontsize,'FontWeight','bold','FontName','Calibri','XTick',[-90 -45 0 45 90],'XTickLabel',...
    {'-90','-45','0','45','90'},'YAxisLocation','right','YTick',...
    [-90 -45 0 45 90]);
xlabel('\theta_1 [\circ]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')
ylabel('\theta_2 [\circ]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')


subplot(1,3,2)
pcolor(t1,t2,a_error_1');
colorbar
shading flat
set(gca,'YaxisLocation','right','Box','off','FontSize',frontsize,'FontWeight','bold','FontName','Calibri','XTick',[-90 -45 0 45 90],'XTickLabel',...
    {'-90','-45','0','45','90'},'YAxisLocation','right','YTick',...
    [-90 -45 0 45 90]);
xlabel('\theta_1 [\circ]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')
ylabel('\theta_2 [\circ]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')

subplot(1,3,3)
pcolor(t1,t2,b_error_1');
colorbar
shading flat
set(gca,'YaxisLocation','right','Box','off','FontSize',frontsize,'FontWeight','bold','FontName','Calibri','XTick',[-90 -45 0 45 90],'XTickLabel',...
    {'-90','-45','0','45','90'},'YAxisLocation','right','YTick',...
    [-90 -45 0 45 90]);
xlabel('\theta_1 [\circ]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')
ylabel('\theta_2 [\circ]','FontSize',frontsize,'FontWeight','bold','FontName','Calibri')

%%
% ratio1_1 = sum(sum(theta_error_1<1))/(length(t1)*length(t2))
% ratio1_2 = sum(sum(theta_error_1<5))/(length(t1)*length(t2))
% ratio2_1 = sum(sum(a_error_1<0.02))/(length(t1)*length(t2))
% ratio2_2 = sum(sum(a_error_1<0.1))/(length(t1)*length(t2))
% ratio3_1 = sum(sum(b_error_1<0.02))/(length(t1)*length(t2))
% ratio3_2 = sum(sum(b_error_1<0.1))/(length(t1)*length(t2))