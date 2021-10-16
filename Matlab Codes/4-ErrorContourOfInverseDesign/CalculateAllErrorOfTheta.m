clc
clear all
close all

% Procedures
% 1. Design space is [theta1 theta2] from -90°to 90°
% 2. Need to input a,b,thickness ratio, temperature
% 3. Import all the target shapes within the design space into the inverse design process and get required parameters and cauculate the error 
% 4. Save theta_error_1,a_error_1,b_error_1 into '.mat' format
% 5. Plot the contour using 'PLotContourOfThetaError.m' 
% It takes about 2hours to run this code.

% -------------- input parameter ---------------------------------------
% --- Geometry ---
t_total = 0.002;            % total thickness
t1_ratio = 1/2;             % thickness ratio of the first layer

a = 0.030;                  % initial length [m]
b = 0.010;                  % initial width  [m]
Temperature = 85;

% ----------------------------------------------------------------------
ii = 0;
jj = 0;
for theta1_step = -90:1:90
    ii = ii+1;
    jj = 0;
    for theta2_step = -90:1:90
        theta = [theta1_step  theta2_step]*pi/180;
        [xx,yy,zz,Pcor,Dtype] = DeformedShape(theta,a,b, t1_ratio,t_total, Temperature);
        if Dtype==0
            [theta1_pred,theta2_pred,a_pred,b_pred] = InPlaneInverseDesign(xx,yy,zz,Pcor,Temperature);
        else
            [theta1_pred,theta2_pred,a_pred,b_pred] = OutOfPlaneInverseDesign(xx,yy,zz,Pcor,t1_ratio,Temperature);
        end
        [theta_e,a_e,b_e,num_pred] = CalculateError(theta1_pred,theta2_pred,a_pred,b_pred,theta,a,b,t1_ratio,Temperature);
        jj = jj+1;
        theta_error_1(ii,jj) = theta_e;
        a_error_1(ii,jj) = a_e;
        b_error_1(ii,jj) = b_e;
        [theta1_step    theta2_step]
    end
end

