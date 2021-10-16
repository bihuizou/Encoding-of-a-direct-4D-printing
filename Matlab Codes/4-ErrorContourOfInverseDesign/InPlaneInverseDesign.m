function [theta1_pred,theta2_pred,a_pred,b_pred] = InPlaneInverseDesign(xx,yy,zz,Pcor,Temperature)
%INPLANEINVERSEDESIGN
% loading map for in_plane deformation at Temperature
load(['inplane_map_',num2str(Temperature),'.mat']);
theta1_pred = [];
theta2_pred = [];
a_pred = [];
b_pred = [];
%    2 -----------  3
%       \          \ 
%        \          \ 
%       1  ----------- 4 

% find theta
dt = 1; % step of theta
j = 0;
for theta_i = -90:dt:90
    [~,i] = min(abs(Angle-theta_i));  % find the location
    % read strains from the map
    ex_r = ex_map(i);
    ey_r = ey_map(i);
    exy_r = exy_map(i);
    % calculate the trial values of a,b
    a_trial = Pcor(2,1)/(ex_r+1);
    b_trial1 = Pcor(4,2)/(ey_r+1);
    b_trial2 = Pcor(4,1)/exy_r;
    if (abs(exy_r)<1e-10) && (abs(Pcor(4,1))<1e-10)
        a_pred(j+1) = a_trial;
        b_pred(j+1) = b_trial1;
        theta1_pred(j+1) = theta_i;
        theta2_pred(j+1) = theta_i;
        j = j+1;
    else
        % check if admissible
        if abs(b_trial1-b_trial2)<1e-10
            % a_p, b_p, t_p are predicted values for a,b,theta
            a_pred(j+1) = a_trial;
            b_pred(j+1) = b_trial1;
            theta1_pred(j+1) = theta_i;
            theta2_pred(j+1) = theta_i;
            j = j+1;
        end
    end
end

end

