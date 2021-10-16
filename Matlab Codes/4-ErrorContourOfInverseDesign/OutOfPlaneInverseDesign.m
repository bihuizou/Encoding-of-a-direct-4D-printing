function [theta1_pred,theta2_pred,a_pred,b_pred] = OutOfPlaneInverseDesign(xx,yy,zz,Pcor,t1_ratio,Temperature)
%OUTOFPLANEINVERSEDESIGN 

%%%%%%%%%% without plotting %%%%%%%%%%%%%

% loading map
switch t1_ratio
    case 0.5
        load(['outplane_6coef_map_1_2_',num2str(Temperature),'.mat']);
    case 1/3
        load(['outplane_6coef_map_1_3_',num2str(Temperature),'.mat']);
    case 1/4
        load(['outplane_6coef_map_1_4_',num2str(Temperature),'.mat']);
    case 2/3
        load(['outplane_6coef_map_2_3_',num2str(Temperature),'.mat']);
    case 3/4
        load(['outplane_6coef_map_3_4_',num2str(Temperature),'.mat']);      
end


[K,H,P1,P2,D1,D2] = surfcurvature(xx,yy,zz);
% K - Gaussian curvature;  H - Mean curvature;  
% P1 - max curvature at each point; D1 - corresponding direction
% P2 - min curvature at each point; D2 - corresponding direction

P1max = max(max(P1));
P2min = min(min(P2));
R_tube = -1/P2min;
R_hole = 1/P1max+R_tube;
cos_theta = -P1*R_hole./(1+P1*R_tube);
R_latitude = R_hole+R_tube.*cos_theta;


xx(isnan(P2)) = NaN;
yy(isnan(P2)) = NaN;
zz(isnan(P2)) = NaN;

% calculate  principal curvatures of 4 corner points - pk1(i),pk2(i)
for i = 1:4
    px(i) = Pcor(i,1);
    py(i) = Pcor(i,2);
    pz(i) = Pcor(i,3);
    pdist = (px(i)-xx).^2+(py(i)-yy).^2+(pz(i)-zz).^2;
    [pxnum(i) pynum(i)] = find(pdist==min(min(pdist)));   
    pk1(i) = abs(P1(pxnum(i),pynum(i)));
    pk2(i) = abs(P2(pxnum(i),pynum(i)));
    
    pr1(i) = abs(R_latitude(pxnum(i),pynum(i)));
end


pk1(5) = abs(P1max);
pk2(5) = abs(P2min);


% define M
zm = pz(2);  % zM = zB

% judge the type
% type1- 1: A,B at different sides 
%        2 :A,B at the same side 
if pz(1)*pz(2)<0
    type1 = 1;
else
    type1 = 2;
end;

% type2- 1: A is higher than B (fai>0)
%        2: B is higher than A (fai<0)
%        3: A and B are at the same height (fai=0)

if pz(1)>pz(2)
    type2 = 1; % fai_pred >0
else
    if pz(1)<pz(2)
        type2 = 2; % fai_pred < 0
    else
        type2 = 3; % fai_pred = 0
    end
end

% flip the point
if type1==1
    pz(2) = -pz(2); % B'=-B
    pz(3) = -pz(3); % C'=-C
else
    pz(3) = -pz(3); % C'=-C
    pz(4) = -pz(4); % D'=-D
end

L_AC = sqrt((px(1)-px(3))^2+(py(1)-py(3))^2);
L_BD = sqrt((px(2)-px(4))^2+(py(2)-py(4))^2);

pk1(1:4) = 1./pr1(1:4);

%% check the central angle is > or < than pi
% find the center of C1 in torus
syms Cx Cy
ra = pr1(1);
rb = pr1(2);
r1 = 1/pk1(5);
sol = solve((Cx-px(1))^2+(Cy-py(1))^2-ra^2==0 , (Cx-px(2))^2+(Cy-py(2))^2-rb^2==0);
O1(:,1) = double(sol.Cx);
O1(:,2) = double(sol.Cy);
Odist = abs((O1(:,1)-px(3)).^2+(O1(:,2)-py(3)).^2-ra^2);
O1 = O1(  find( Odist==min(Odist)),   :);
OA = [px(1)-O1(1),py(1)-O1(2)];
OB = [px(2)-O1(1),py(2)-O1(2)];
OC = [px(3)-O1(1),py(3)-O1(2)];
OD = [px(4)-O1(1),py(4)-O1(2)];
crossCA = OC(1)*OA(2)-OA(1)*OC(2);
crossBD = OB(1)*OD(2)-OD(1)*OB(2);


%% determine theta
m0 = L_AC/2/ra;
n0 = L_BD/2/rb;
if abs(m0)>1 
    m0 = fix(m0);
end
if abs(n0)>1
    n0 = fix(n0);
end

if crossCA>0
    t1 = 2*asin(m0); % theta1
else 
    t1 = 2*pi-2*asin(m0); % theta1
end
if crossBD>0
    t2 = 2*asin(n0); % theta2
else
    t2 = 2*pi-2*asin(n0); % theta2
end
t3 = (t1+t2)/2;         % theta3



% arc MB
MN = t1/pk1(2);
arcMB = t3*rb;
MB = arcMB*(abs(pk1(2))/P1max);

% calculate central angle of view--z-xy
r2 = 1/pk2(5);
%
if ra>(r1+r2)
    t4 = pi-asin(abs(pz(1))/r2); % theta4
else 
    t4 = asin(abs(pz(1))/r2); % theta4
end
if rb>(r1+r2)
    t5 = pi-asin(abs(zm)/r2);    % theta5
else
    t5 = asin(abs(zm)/r2);    % theta5
end

if type1==1
    t6 = t4+t5;
else
    t6 = abs(t4-t5);
end
arcMA = r2*t6;
MA = arcMA;


if type2==1
    fai_pred = atan(arcMA/MB);
else
    fai_pred = atan(-arcMA/MB);
end



%%% find theta1, theta2, a, b
% calculate kx,ky,kxy from k1,k2,fai_pred
p = (P1max+P2min)/2;
rr = (P1max-P2min)/2;
kx_p = p + rr*cos(2*fai_pred);
ky_p = p - rr*cos(2*fai_pred);
kxy_p = 2*rr*sin(-2*fai_pred);


%% find the corresponding theta1_p, theta2_p from the map
% 1. first possibility
dist = sqrt( (kx_map-kx_p).^2 + (ky_map-ky_p).^2 + (kxy_map-kxy_p).^2  );
[index_x index_y] = find(dist==min(min(dist))); 
dist_min = dist(index_x(1),index_y(1));
[all_posx all_posy] = find(abs(dist-dist_min)<0.1);  % find all the possible theta
possible_theta1(:,1) = Angle(all_posx);
possible_theta1(:,2) = Angle(all_posy);
% 2. second possibility
dist = sqrt( (kx_map-(-kx_p)).^2 + (ky_map-(-ky_p)).^2 + (kxy_map-(-kxy_p)).^2  );
[index_x index_y] = find(dist==min(min(dist))); 
dist_min = dist(index_x(1),index_y(1));
[all_posx all_posy] = find(abs(dist-dist_min)<0.1);  % find all the possible theta
possible_theta2(:,1) = Angle(all_posx);
possible_theta2(:,2) = Angle(all_posy);
possible_theta = [possible_theta1;possible_theta2];

% by this, we can get two groups of (theta1,theta2) that satisfies the target curvatures,

% So we need to analyze the in-plane case to determine exact(theta1,theta2) , and the initial length a and b.
% 1. calculate 3 strains
% calculate the length of 2D flatten plot

if (type2==1)   % zA>zB
    if t1>t2
        Lx = t1/P1max;
        Ly = 2*t5/abs(P2min);
    else
        Lx = t2/P1max;
        Ly = 2*t5/abs(P2min);
    end
else          % zA<zB
    if t1<t2
        Lx = t2/P1max;
        Ly = 2*t4/abs(P2min);
    else
        Lx = t1/P1max;
        Ly = 2*t4/abs(P2min);
    end
end


% calculate the coordinates of points ABCD in 2D flatten plate
% Pos_pc --- position in principal curvature coordinate system
if type2==1   % zA>zB
    if t1>t2
        Pos_pc(1,1) = 0;      Pos_pc(1,2) = 0;         % point A
        Pos_pc(2,1) = MB;     Pos_pc(2,2) = -MA;       % point B
        Pos_pc(3,1) = Lx;     Pos_pc(3,2) = Ly-2*MA;   % point C
        Pos_pc(4,1) = Lx-MB;  Pos_pc(4,2) = Ly-MA;     % point D
    else
        Pos_pc(1,1) = 0;        Pos_pc(1,2) = 0;         % point A
        Pos_pc(2,1) = MB;       Pos_pc(2,2) = -MA;       % point B
        Pos_pc(3,1) = 2*MB-Lx;  Pos_pc(3,2) = Ly-2*MA;   % point C
        Pos_pc(4,1) = MB-Lx;    Pos_pc(4,2) = Ly-MA;     % point D
    end
else          % zA<zB
    if type2==2
        if t1<t2
            Pos_pc(1,1) = 0;        Pos_pc(1,2) = 0;        % point A
            Pos_pc(2,1) = MB;       Pos_pc(2,2) = MA;       % point B
            Pos_pc(3,1) = 2*MB-Lx;  Pos_pc(3,2) = Ly;       % point C
            Pos_pc(4,1) = MB-Lx;    Pos_pc(4,2) = Ly-MA;    % point D
        else
            Pos_pc(1,1) = 0;        Pos_pc(1,2) = 0;        % point A
            Pos_pc(2,1) = MB;       Pos_pc(2,2) = MA;       % point B
            Pos_pc(3,1) = Lx;       Pos_pc(3,2) = Ly;       % point C
            Pos_pc(4,1) = Lx-MB;    Pos_pc(4,2) = Ly-MA;    % point D
        end
    else      % zA = zB
        Pos_pc(1,1) = 0;        Pos_pc(1,2) = 0;        % point A
        Pos_pc(2,1) = MB;       Pos_pc(2,2) = 0;        % point B
        Pos_pc(3,1) = MB;       Pos_pc(3,2) = Ly;       % point C
        Pos_pc(4,1) = 0;        Pos_pc(4,2) = Ly;       % point D
    end
end

% transform coordinatse into structural coordinate system

for i = 1:4
    Pos_s(i,1) = cos(fai_pred)*Pos_pc(i,1)-sin(fai_pred)*Pos_pc(i,2);
    Pos_s(i,2) = sin(fai_pred)*Pos_pc(i,1)+cos(fai_pred)*Pos_pc(i,2);
end


 
% check whethe the possible_theta is admissible
j = 0;

preset_error = 1e-2;
JudgeWhetherNaN = 0;
while (JudgeWhetherNaN == 0) &&  (preset_error<=10)
    for i = 1:size(possible_theta,1)
        % for i = 1:1
        [~,index_x] = min(abs(Angle-possible_theta(i,1)));
        [~,index_y] = min(abs(Angle-possible_theta(i,2)));
        % read strains from the map
        ex_r = ex_map(index_x,index_y);
        ey_r = ey_map(index_x,index_y);
        exy_r = exy_map(index_x,index_y);
        
        % calculate the trial values of a,b
        a_trial = Pos_s(2,1)/(ex_r+1);
        b_trial1 = Pos_s(4,2)/(ey_r+1);
        b_trial2 = Pos_s(4,1)/exy_r;
        
        % check if admissible
        if (abs(exy_r)<1e-5) && (abs(Pos_s(4,1))<1e-3)    % without shear
            a_pred(j+1) = a_trial;
            b_pred(j+1) = b_trial1;
            theta1_pred(j+1) = possible_theta(i,1);
            theta2_pred(j+1) = possible_theta(i,2);
            JudgeWhetherNaN = 1;
            j = j+1;
        else
            %         if abs(Pos_s(4,2)*exy_r-Pos_s(4,1)*(ey_r+1))<1e-5
            if abs(b_trial1-b_trial2)<preset_error
                %      if abs(exy_trial-exy_r)<1e-3
                % a_p, b_p, t_p are predicted values for a,b,theta
                a_pred(j+1) = a_trial;
                b_pred(j+1) = b_trial1;
                theta1_pred(j+1) = possible_theta(i,1);
                theta2_pred(j+1) = possible_theta(i,2);
                JudgeWhetherNaN = 1;
                j = j+1;
            end
        end
    end
    preset_error = preset_error*10;
    if (preset_error==10 ) && (abs(exy_r)<1e-5)
                    a_pred(j+1) = a_trial;
            b_pred(j+1) = b_trial1;
            theta1_pred(j+1) = possible_theta(i,1);
            theta2_pred(j+1) = possible_theta(i,2);
            JudgeWhetherNaN = 1;
            j = j+1;
    end
end

end

