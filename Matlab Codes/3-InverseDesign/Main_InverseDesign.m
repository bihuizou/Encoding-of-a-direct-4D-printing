clc
clear all
close all

% Procedures
% 1. Input [geometrical parameters, thickness ratio, activation temperature] and generate the target shape 
% 2. Import the target shape into the inverse design process and get the 

% -------------- input parameter ---------------------------------------
% --- Geometry ---
t_total = 0.002;            % total thickness
t1_ratio = 2/3;             % thickness ratio of the first layer
theta = [0 90]*pi/180;      % angle of each layer
a = 0.030;                  % initial length [m]
b = 0.010;                  % initial width  [m]
Temperature = 85;

% ----------------------------------------------------------------------

%% generate the coordinates of target shape
[xx,yy,zz,Pcor,Dtype] = DeformedShape(theta,a,b, t1_ratio,t_total, Temperature);

%% inverse design 
if Dtype==0
    [theta1_pred,theta2_pred,a_pred,b_pred] = InPlaneInverseDesign(xx,yy,zz,Pcor,Temperature);
else
    [theta1_pred,theta2_pred,a_pred,b_pred] = OutOfPlaneInverseDesign(xx,yy,zz,Pcor,t1_ratio,Temperature);
end

%% calculate the error between input values and predicted values
[theta_error,a_error,b_error,num_pred] = CalculateError(theta1_pred,theta2_pred,a_pred,b_pred,theta,a,b,t1_ratio,Temperature);

%% output
disp(['The predicted theta_1 is ',num2str(theta1_pred(num_pred)),'°'])
disp(['The predicted theta_2 is ',num2str(theta2_pred(num_pred)),'°'])
disp(['The predicted length a is ',num2str(a_pred(num_pred)),'m'])
disp(['The predicted length b is ',num2str(b_pred(num_pred)),'m'])