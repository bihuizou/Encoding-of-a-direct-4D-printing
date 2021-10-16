clc
clear all
close all

% Procedures
% 1. Design space is [a b] from 10mm to 40mm
% 2. Need to input theta1, theta2, thickness ratio, temperature
% 3. Import all the target shapes within the design space into the inverse design process and get required parameters and cauculate the error 
% 4. Save theta_error_1,a_error_1,b_error_1 into '.mat' format
% 5. Plot the contour using 'PLotContourOfLengthError.m' 


% -------------- input parameter ---------------------------------------
% --- Geometry ---
t_total = 0.002;            % total thickness
t1_ratio = 2/3;             % thickness ratio of the first layer
theta = [0 45]*pi/180;
Temperature = 85;

% ----------------------------------------------------------------------
ii = 0;
jj = 0;
a = 0.01:0.001:0.04;
b = 0.01:0.001:0.04;
theta_error_1 = zeros(length(a),length(b))*NaN;
a_error_1 = zeros(length(a),length(b))*NaN;
b_error_1 = zeros(length(a),length(b))*NaN;

for i = 1 :length(a)
    ii = ii+1;
    jj = 0;
    for j = 1 :i

        [xx,yy,zz,Pcor,Dtype] = DeformedShape(theta,a(i),b(j), t1_ratio,t_total, Temperature);
        if Dtype==0
            [theta1_pred,theta2_pred,a_pred,b_pred] = InPlaneInverseDesign(xx,yy,zz,Pcor,Temperature);
        else
            [theta1_pred,theta2_pred,a_pred,b_pred] = OutOfPlaneInverseDesign(xx,yy,zz,Pcor,t1_ratio,Temperature);
        end
        [theta_e,a_e,b_e,num_pred] = CalculateError(theta1_pred,theta2_pred,a_pred,b_pred,theta,a(i),b(j),t1_ratio,Temperature);
        jj = jj+1;
        theta_error_1(ii,jj) = theta_e;
        a_error_1(ii,jj) = a_e;
        b_error_1(ii,jj) = b_e;
        [a(i)  b(j)]
    end
end

