function [xx,yy,zz,Pcor,Dtype] = DeformedShape(theta,a,b, t1_ratio,t_total, Temperature)
%DEFORMEDSHAPE

%%%%%%%%%%%  without plotting   %%%%%%%%%%%%%%%%%


% --- Material ---
% E1,E2,G12: modulus [Pa];   v12,v21: Poisson's Ratio;    e12_T:  thermal strain 
load(['MaterialData_',num2str(Temperature),'.mat']);

% --- Geometry ---
thick1 = t_total*t1_ratio;
thick2 = t_total-thick1;
t = [thick1 thick2];               % thickness of each layer
n = length(t);                     % number of layers


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
if abs(k(3))<1e-6 
    k(3) = 0;
end

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
            if term<0 
                judge = 0;
            end
        end
        if judge ==0 
            zz(i,j) = NaN;
        end
    end
end




%% judge whether the deformation is in-plane or out-of-plane
if (abs(k(1))<1e-6) &  (abs(k(2))<1e-6) & (abs(k(3))<1e-6) 
    Dtype = 0;  % inplane
else
    Dtype = 1;  % out-of-plane
end

%% obtain the coornidate of 4 corner points of target shape
if Dtype==0
    % Pos_t - coordinates of target shape
    Pcor(1,1) = rect(1,1); Pcor(1,2) = rect(1,2); Pcor(1,3) = 0;% point A
    Pcor(2,1) = rect(4,1); Pcor(2,2) = rect(4,2); Pcor(2,3) = 0 % point B
    Pcor(3,1) = rect(3,1); Pcor(3,2) = rect(3,2); Pcor(3,3) = 0 % point C
    Pcor(4,1) = rect(2,1); Pcor(4,2) = rect(2,2); Pcor(4,3) = 0 % point D
else
    Pcor_2d(3,1) = x2(1,1);                  Pcor_2d(3,2) = y2(1,1);                  %C
    Pcor_2d(4,1) = x2(length(x),1);          Pcor_2d(4,2) = y2(length(x),1);          %D
    Pcor_2d(1,1) = x2(length(x),length(y));  Pcor_2d(1,2) = y2(length(x),length(y));  %A
    Pcor_2d(2,1) = x2(1,length(y));          Pcor_2d(2,2) = y2(1,length(y));          %B
    for i = 1:4
        v = theta1.*(Pcor_2d(i,1)-x_l)/L1;
        u = (theta2.*(Pcor_2d(i,2)-y_l)/L2)+pi-theta2/2;
        Pcor(i,1) = (R+r*cos(u)).*cos(v);
        Pcor(i,2) = (R+r*cos(u)).*sin(v);
        Pcor(i,3) = r*sin(u)+0*v;
    end
end
end

