clc
clear all
close all

% -------------- input parameter ---------------------------------------
% --- Geometry ---
t1_ratio = 1/2;             % thickness ratio of the first layer

% --- Material ---
% E1,E2,G12: modulus [Pa];   v12,v21: Poisson's Ratio;    e12_T:  thermal strain 
% E1=Mat(1); v12=Mat(2); v21=Mat(3); E2=Mat(4); G12=Mat(5);

% Mat = [8.26e6 0.5 0.4 6.61e6 0.79e6];    e12_T = [-0.05787; 0.00681;  0];  % T=65; % T=65;
% Mat = [5.82e6 0.5 0.4 4.66e6 0.78e6];  e12_T = [-0.10680; 0.02568; 0];   % T=75;% T=75;
Mat = [4.87e6 0.5 0.4 3.90e6 0.67e6];  e12_T = [-0.14539; 0.02997; 0];  % T=85;% T=85;

% ----------------------------------------------------------------------


%%
t_total = 0.002;           % total thickness
t1 = t_total*t1_ratio;
t2 = t_total-t1;

angle1 = -90:1:90;
angle2 = -90:1:90;



for i = 1:length(angle1)
    for j = 1:length(angle2)
        [k1(i,j),k2(i,j),fai(i,j),ex(i,j),ey(i,j),exy(i,j),kx(i,j),ky(i,j),kxy(i,j),C(i,j)] = cal_k(angle1(i),angle2(j),t1,t2,e12_T,Mat);
        fai(i,j) = fai(i,j)*180/pi;
    end
end

%% find the modes
% in-plane  kx=ky=kxy=0
% mode1: extension   kx=ky=kxy=0, exy=0
% mode2:  shearing   kx=ky=kxy=0 
% out-of-plane 
% mode3: bending     kxy=0
% mode4: twisting    otherwise

[inplane_extend_x,inplane_extend_y] = find( (abs(kxy)<1e-5) & (abs(kx)<1e-5) & (abs(ky)<1e-5) & (abs(exy)<1e-5) );
[inplane_shear_x,inplane_shear_y] = find( (abs(kxy)<1e-5) & (abs(kx)<1e-5) & (abs(ky)<1e-5) & (abs(exy)>1e-5) );
[outplane_bend_x,outplane_bend_y] = find( (abs(kxy)<1e-5) );
[outplane_twist_x,outplane_twist_y] = find( (abs(kxy)>1e-5) );
inplane_extend_x(inplane_extend_x==181) = 180;
inplane_extend_y(inplane_extend_y==181) = 180;





%%
figure(1)
subplot(2,3,1)
pcolor(angle1,angle2,ex')
xlabel('\theta_1','FontWeight','bold','FontName','Calibri');
ylabel('\theta_2','FontWeight','bold','FontName','Calibri');
colorbar
shading flat
hold on
contour(angle1,angle2,kxy',[0 0],'--','color','k','LineWidth',1.2)
plot(angle1(inplane_shear_x),angle2(inplane_shear_y),'color',[204/255,0/255,255/255],'LineWidth',1.5)
plot(angle1(inplane_extend_x),angle2(inplane_extend_y),'s','MarkerFaceColor',[192/255,191/255,191/255],'MarkerEdgeColor','k','markersize',10)
set(gca,'FontName','Calibri','FontWeight','bold','XTick',...
    [-90 -45 0 45 90],'YTick',[-90 -45 0 45 90]);
title('e0_x')


subplot(2,3,2)
pcolor(angle1,angle2,ey')
xlabel('\theta_1','FontWeight','bold','FontName','Calibri');
ylabel('\theta_2','FontWeight','bold','FontName','Calibri');
colorbar
shading flat
hold on
contour(angle1,angle2,kxy',[0 0],'--','color','k','LineWidth',1.2)
plot(angle1(inplane_shear_x),angle2(inplane_shear_y),'color',[204/255,0/255,255/255],'LineWidth',1.5)
plot(angle1(inplane_extend_x),angle2(inplane_extend_y),'s','MarkerFaceColor',[192/255,191/255,191/255],'MarkerEdgeColor','k','markersize',10)
set(gca,'FontName','Calibri','FontWeight','bold','XTick',...
    [-90 -45 0 45 90],'YTick',[-90 -45 0 45 90]);
title('e0_y')

subplot(2,3,3)
pcolor(angle1,angle2,exy')
xlabel('\theta_1','FontWeight','bold','FontName','Calibri');
ylabel('\theta_2','FontWeight','bold','FontName','Calibri');
colorbar
shading flat
hold on
contour(angle1,angle2,kxy',[0 0],'--','color','k','LineWidth',1.2)
plot(angle1(inplane_shear_x),angle2(inplane_shear_y),'color',[204/255,0/255,255/255],'LineWidth',1.5)
plot(angle1(inplane_extend_x),angle2(inplane_extend_y),'s','MarkerFaceColor',[192/255,191/255,191/255],'MarkerEdgeColor','k','markersize',10)
set(gca,'FontName','Calibri','FontWeight','bold','XTick',...
    [-90 -45 0 45 90],'YTick',[-90 -45 0 45 90]);
title('e0_x_y')

subplot(2,3,4)
pcolor(angle1,angle2,kx')
xlabel('\theta_1','FontWeight','bold','FontName','Calibri');
ylabel('\theta_2','FontWeight','bold','FontName','Calibri');
colorbar
shading flat
hold on
contour(angle1,angle2,kxy',[0 0],'--','color','k','LineWidth',1.2)
plot(angle1(inplane_shear_x),angle2(inplane_shear_y),'color',[204/255,0/255,255/255],'LineWidth',1.5)
plot(angle1(inplane_extend_x),angle2(inplane_extend_y),'s','MarkerFaceColor',[192/255,191/255,191/255],'MarkerEdgeColor','k','markersize',10)
set(gca,'FontName','Calibri','FontWeight','bold','XTick',...
    [-90 -45 0 45 90],'YTick',[-90 -45 0 45 90]);
title('k_x')

subplot(2,3,5)
pcolor(angle1,angle2,ky')
xlabel('\theta_1','FontWeight','bold','FontName','Calibri');
ylabel('\theta_2','FontWeight','bold','FontName','Calibri');
colorbar
shading flat
hold on
contour(angle1,angle2,kxy',[0 0],'--','color','k','LineWidth',1.2)
plot(angle1(inplane_shear_x),angle2(inplane_shear_y),'color',[204/255,0/255,255/255],'LineWidth',1.5)
plot(angle1(inplane_extend_x),angle2(inplane_extend_y),'s','MarkerFaceColor',[192/255,191/255,191/255],'MarkerEdgeColor','k','markersize',10)
set(gca,'FontName','Calibri','FontWeight','bold','XTick',...
    [-90 -45 0 45 90],'YTick',[-90 -45 0 45 90]);
title('k_y')

subplot(2,3,6)
pcolor(angle1,angle2,kxy')
xlabel('\theta_1','FontWeight','bold','FontName','Calibri');
ylabel('\theta_2','FontWeight','bold','FontName','Calibri');
colorbar
shading flat
hold on
contour(angle1,angle2,kxy',[0 0],'--','color','k','LineWidth',1.2)
plot(angle1(inplane_shear_x),angle2(inplane_shear_y),'color',[204/255,0/255,255/255],'LineWidth',1.5)
plot(angle1(inplane_extend_x),angle2(inplane_extend_y),'s','MarkerFaceColor',[192/255,191/255,191/255],'MarkerEdgeColor','k','markersize',10)
set(gca,'FontName','Calibri','FontWeight','bold','XTick',...
    [-90 -45 0 45 90],'YTick',[-90 -45 0 45 90]);
title('k_x_y')
