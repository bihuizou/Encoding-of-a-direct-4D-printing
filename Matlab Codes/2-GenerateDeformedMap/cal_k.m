
function [k1,k2,fai,ex,ey,exy,kx,ky,kxy,C] = cal_k(angle1,angle2,t1,t2,e12_T,Mat)

% Properties
n = 2;                % number of layers
t = [t1 t2];          % thickness of each layer
theta = [angle1 angle2]*pi/180;     % angle of each layer
% E1=Mat(1); v12=Mat(2); v21=Mat(3); E2=E1*v21/v21; G12=Mat(5);
E1=Mat(1); v12=Mat(2); v21=Mat(3); E2=Mat(4); G12=Mat(5);


Q = [ E1/(1-v12*v21)     E1*v21/(1-v12*v21)       0;
     E1*v21/(1-v12*v21)    E2/(1-v12*v21)         0;
         0                       0               G12  ];



     
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
    T = [cos(theta(i))^2  sin(theta(i))^2  2*sin(theta(i))*cos(theta(i));
         sin(theta(i))^2  cos(theta(i))^2  -2*sin(theta(i))*cos(theta(i));
         -sin(theta(i))*cos(theta(i))  sin(theta(i))*cos(theta(i)) cos(theta(i))^2-sin(theta(i))^2;];
    Qbar = inv(T)*Q*(inv(T))';
    
    exy_T = (T')*e12_T;
    
    N_T = N_T+Qbar*exy_T*(z(i+1)-z(i));
    M_T = M_T+1/2*Qbar*exy_T*(z(i+1)^2-z(i)^2);
    
    A = A+Qbar*(z(i+1)-z(i));
    B = B+1/2*Qbar*(z(i+1)^2-z(i)^2);
    D = D+1/3*Qbar*(z(i+1)^3-z(i)^3);
    
end

f = [A B;B D];
sol = inv([A B;B D])*[N_T;M_T];
e0 = sol(1:3,1);
k  = sol(4:6,1);
% e0 = sol(1);
% k  = sol(4);

% principle curvature
C = 1/2*(k(1)+k(2));
R = sqrt(  ((k(1)-k(2))/2)^2  +  (k(3)/2)^2      );
fai = -1/2*atan( k(3)/(k(1)-k(2))   );
k1 = C+R;
k2 = C-R;
kx = k(1);
ky = k(2);
kxy = k(3);
ex = e0(1);
ey = e0(2);
exy = e0(3);



end

