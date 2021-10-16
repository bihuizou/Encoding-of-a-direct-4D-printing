clc
clear all
close all


% -------------- input parameter ---------------------------------------
% --- Geometry ---
t1_ratio = 1/2;             % thickness ratio of the first layer
theta = [62 -64]*pi/180;      % angle of each layer
a = 0.030;                  % initial length [m]
b = 0.010;                  % initial width  [m]

% --- Material ---
% E1,E2,G12: modulus [Pa];   v12,v21: Poisson's Ratio;    e12_T:  thermal strain 

% E1 = 8.26e6; v12 = 0.5; v21 = 0.4; E2 = E1*v21/v12; G12 = 0.79e6; e12_T = [-0.05787; 0.00681; 0]; % T=65
% E1 = 5.82e6; v12 = 0.5; v21 = 0.4; E2 = E1*v21/v12; G12 = 0.78e6; e12_T = [-0.10680; 0.02568; 0]; % T=75
E1 = 4.87e6; v12 = 0.5; v21 = 0.4; E2 = E1*v21/v12; G12 = 0.67e6; e12_T = [-0.14539; 0.02997; 0]; % T=85  

% ----------------------------------------------------------------------

t_total = 0.002;            % total thickness
t1 = t_total*t1_ratio;
t2 = t_total-t1;
t = [t1 t2];               % thickness of each layer
n = length(t);              % number of layers

%% Classical laminate theory to calculate the strain and curvature

Q = [ E1/(1-v12*v21)       E1*v21/(1-v12*v21)       0;
     E1*v21/(1-v12*v21)     E2/(1-v12*v21)          0;
         0                        0                G12];
     
% calculate the z-coordanates 
m = sum(t)/2;
M = zeros(1,n+1);
for i=1:n
    M(i+1) = M(i)+t(i);
end
z = M-m;     

% thermal resultant force and moment
N_T = 0;
M_T = 0;

A = 0;
B = 0;
D = 0;

% assemble 
for i=1:n
    % transformation matrix
    T(:,:,i) = [cos(theta(i))^2  sin(theta(i))^2  2*sin(theta(i))*cos(theta(i));
         sin(theta(i))^2  cos(theta(i))^2  -2*sin(theta(i))*cos(theta(i));
         -sin(theta(i))*cos(theta(i))  sin(theta(i))*cos(theta(i)) cos(theta(i))^2-sin(theta(i))^2;];
    Qbar(:,:,i) = inv(T(:,:,i))*Q*(inv(T(:,:,i)))';
    
    exy_T = (T(:,:,i)')*e12_T;
    
    N_T = N_T+Qbar(:,:,i)*exy_T*(z(i+1)-z(i));
    M_T = M_T+1/2*Qbar(:,:,i)*exy_T*(z(i+1)^2-z(i)^2);
    
    A = A+Qbar(:,:,i)*(z(i+1)-z(i));
    B = B+1/2*Qbar(:,:,i)*(z(i+1)^2-z(i)^2);
    D = D+1/3*Qbar(:,:,i)*(z(i+1)^3-z(i)^3);
    
end

ABD = [A B;B D];
invABD = inv([A B;B D]);
sol = invABD*[N_T;M_T];
e0 = sol(1:3,1);            % strain at mid-plane
k  = sol(4:6,1);            % curvature

%% Use principal curvature to visualize the deformation
% principle curvature
Center = 1/2*(k(1)+k(2));
R = sqrt(  ((k(1)-k(2))/2)^2  +  (k(3)/2)^2    );
fai = 1/2*atan( k(3)/(k(1)-k(2))   );
k1 = Center+R;
k2 = Center-R;

%%%% out of plane deflection
% plot the in-plane deformation
% 1. initial shape (x,y) in structural coordinate system
step = b/20;
x = 0:step:a;
y = 0:step:b;
x = x';

u0 = e0(1)*x+e0(3)*y;     % displacement
v0 = e0(2)*y+0*x;         % displacement
x1 = u0+x;
y1 = v0+y;
% (x1,y1) deformed shape in structural coordinate system

if (k(1)==0 & k(2)==0 & k(3)==0)    % in-plane deformation
    x2 = x1;
    y2 = y1;
else                                % out-of-plane deformation
    % transformation to principle curvature direction
    fai = -fai;
    for i=1:length(x)
        for j=1:length(y)
            x2(i,j) = x1(i,j).*cos(fai)+y1(i,j).*sin(fai);
            y2(i,j) = -x1(i,j).*sin(fai)+y1(i,j).*cos(fai);
        end
    end
end
% (x2,y2) deformed shape in principal curvature's coordinate system

% new plate
x_l = min(min(x2(:,:)));
x_r = max(max(x2(:,:)));
y_l = min(min(y2(:,:)));
y_r = max(max(y2(:,:)));

L1 = x_r - x_l;
L2 = y_r - y_l;

x3 = x_l:L1/500:x_r;
x3 = x3';
y3 = y_l:L2/500:y_r;
% (x3,y3) imaginary plane in structural coordinate system


theta1 = k1*L1;
theta2 = -k2*L2;


r = -1/k2;
R = 1/k1;
R = R+r;

v = theta1.*(x3-x_l)/L1;
u = (theta2.*(y3-y_l)/L2)+pi-theta2/2;
[u,v] = meshgrid(u,v);

xx0 = (R+r*cos(u)).*cos(v);
yy0 = (R+r*cos(u)).*sin(v);
zz0 = r*sin(u)+0*v;


%% subtract the void points from the imaginary plane to get the actural shape
xx = xx0;
yy = yy0;
zz = zz0;

%% judge whether the point is in the plate
x0 = (x_r+x_l)/2;
y0 = (y_r+y_l)/2;  % point inside the plate

% line
% 4 points at the corner
rect(1,1) = x2(1,1); rect(1,2) = y2(1,1);
rect(2,1) = x2(1,end); rect(2,2) = y2(1,end);
rect(3,1) = x2(end,end); rect(3,2) = y2(end,end);
rect(4,1) = x2(end,1); rect(4,2) = y2(end,1);
rect(5,:) = rect(1,:);

% coefficient of line equations  Ax+By+C=0  
%A-coef(:,1);B-coef(:,2);C-coef(:,3);
for i=1:4
    coef(i,1) = rect(i+1,2)-rect(i,2);
    coef(i,2) = rect(i,1)-rect(i+1,1);
    coef(i,3) = rect(i+1,1)*rect(i,2)-rect(i,1)*rect(i+1,2);
    term1(i) = coef(i,1)*x0+coef(i,2)*y0+coef(i,3);
end


for i =1:length(x3)
    for j =1:length(y3)
        judge = 1;
        for kk = 1:4
            term = term1(kk)*(  coef(kk,1)*x3(i)+coef(kk,2)*y3(j)+coef(kk,3)    );
            Term(i,j)=term;
            if term<0 
                judge = 0;
            end
        end
        if judge ==0 
            zz(i,j) = NaN;
        end
    end
end
%%
if (k(1)==0 & k(2)==0 & k(3)==0)    % in-plane deformation
    figure(1)
    hold on
    patch(rect(1:4,1),rect(1:4,2),[0;1;1;0],'EdgeColor','none','FaceAlpha',0.4)
    plot([rect(1:4,1);rect(1,1)],[rect(1:4,2);rect(1,2)],'k')
    set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')   
    axis equal
    s.EdgeColor = 'none';
    grid off
    color1 = [0.2,0.1,0.9];
    color2 = [0.1,0.9,0.6];
    colorscale1 = [linspace(color2(1),color1(1),100)',linspace(color2(2),color1(2),100)',linspace(color2(3),color1(3),100)'];
    colorscale2 = [linspace(color1(1),color2(1),100)',linspace(color1(2),color2(2),100)',linspace(color1(3),color2(3),100)'];
    colorscale=[colorscale2 ; colorscale1];
    colormap(colorscale);
else
    figure(1)  % full scale
    surf(xx0,yy0,zz0,'FaceAlpha',0.5)
    hold on  
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    shading interp
    
    figure(2)  % after cutting
    s = surf(xx,yy,zz,'FaceAlpha',0.4)
    axis equal
    shading interp
    hold on
    title('Deformed shape')

    set(gcf,'Color','White');
    set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','White','ycolor','White','zcolor','White')
    gridsmal2=1:25:size(xx,1)/2+1;
    gridsmall=1:25:size(xx,1)/2+1;
 
    s.EdgeColor = 'none';
    grid off
    color1 = [0.2,0.1,0.9];
    color2 = [0.1,0.9,0.6];
    colorscale1 = [linspace(color2(1),color1(1),100)',linspace(color2(2),color1(2),100)',linspace(color2(3),color1(3),100)'];
    colorscale2 = [linspace(color1(1),color2(1),100)',linspace(color1(2),color2(2),100)',linspace(color1(3),color2(3),100)'];
    colorscale=[colorscale2 ; colorscale1];
    colormap(colorscale);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    if sum(isnan(zz(:)))>2000
        a=[];
        b=[];
        % for i=1:size(xx,1)
        for i=1:501
            judge2 = prod(isnan(zz(i,:)));
            if judge2==0
                validIndices=~isnan(zz(i,:));
                validnumbers=zz(i,validIndices);
                colume=zz(i,:);
                number1=find(colume==validnumbers(1));
                number2=find(colume==validnumbers(end));
                a=[a,[i;number1]];
                b=[b,[i;number2]];
            end
        end
        lineX1=[];
        lineY1=[];
        lineZ1=[];
        for i=1:size(a,2)
            % plot3(xx(a(1,i),a(2,i)),yy(a(1,i),a(2,i)),zz(a(1,i),a(2,i)),'r.')
            lineX1=[lineX1,xx(a(1,i),a(2,i))];
            lineY1=[lineY1,yy(a(1,i),a(2,i))];
            lineZ1=[lineZ1,zz(a(1,i),a(2,i))];
        end
        lineX2=[];
        lineY2=[];
        lineZ2=[];
        for i=1:size(b,2)
            % plot3(xx(b(1,i),b(2,i)),yy(b(1,i),b(2,i)),zz(b(1,i),b(2,i)),'r.')
            lineX2=[lineX2,xx(b(1,i),b(2,i))];
            lineY2=[lineY2,yy(b(1,i),b(2,i))];
            lineZ2=[lineZ2,zz(b(1,i),b(2,i))];
        end
        % plot3(lineX1,lineY1,lineZ1,'r','LineWidth',3)
        % plot3(lineX2,lineY2,lineZ2,'r','LineWidth',3)
        
        plot3(lineX1,lineY1,lineZ1,'k')
        plot3(lineX2,lineY2,lineZ2,'k')
        
    else
        plot3(xx(1,:),yy(1,:),zz0(1,:),'k')  % draw the mesh
        plot3(xx(:,1)',yy(:,1)',zz0(:,1)','k')  % draw the mesh
        plot3(xx(end,:),yy(end,:),zz0(end,:),'k')  % draw the mesh
        plot3(xx(:,end)',yy(:,end)',zz0(:,end)','k')  % draw the mesh
    end
    
end



        





