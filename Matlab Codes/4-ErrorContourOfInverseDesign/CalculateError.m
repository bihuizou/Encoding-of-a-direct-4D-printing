function [theta_error,a_error,b_error,num_pred] = CalculateError(theta1_pred,theta2_pred,a_pred,b_pred,theta,a,b,t1_ratio,Temperature)
%CALCULATEERROR

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


%% error of theta
ThetaEqual = 0;
for i = 1:length(theta1_pred)
    if (theta(1)== theta1_pred(i)*pi/180) && (theta(2)== theta2_pred(i)*pi/180)
        num_pred = i;
        ThetaEqual = 1;
        theta_error = 0;
        break;
    else
        [~,index_x] = min(abs(Angle- theta1_pred(i)   ));
        [~,index_y] = min(abs(Angle- theta2_pred(i)   ));
        % read strains from the map
        ex_r(i) = ex_map(index_x,index_y);
        ey_r(i) = ey_map(index_x,index_y);
        exy_r(i) = exy_map(index_x,index_y);
        kx_r(i) = kx_map(index_x,index_y);
        ky_r(i) = ky_map(index_x,index_y);
        kxy_r(i) = kxy_map(index_x,index_y);
    end
end

if ThetaEqual == 0
    [~,index_x] = min(abs(Angle*pi/180- theta(1)   ));
    [~,index_y] = min(abs(Angle*pi/180- theta(2)   ));
    % read strains from the map
    ex_real_r = ex_map(index_x,index_y);
    ey_real_r = ey_map(index_x,index_y);
    exy_real_r = exy_map(index_x,index_y);
    kx_real_r = kx_map(index_x,index_y);
    ky_real_r = ky_map(index_x,index_y);
    kxy_real_r = kxy_map(index_x,index_y);
    
    Sum = 0;
        if abs(ex_real_r)>1e-10
        Sum = ( (ex_r-ex_real_r)./ex_real_r ).^2;
    end
    if abs(ey_real_r)>1e-10
        Sum = Sum + ( (ey_r-ey_real_r)./ey_real_r ).^2;
    end
    if abs(exy_real_r)>1e-10
        Sum = Sum + ( (exy_r-exy_real_r)./exy_real_r ).^2;
    end
    if abs(kx_real_r)>1e-10
        Sum = Sum + ( (kx_r-kx_real_r)./kx_real_r ).^2;
    end
    if abs(ky_real_r)>1e-10
        Sum = Sum + ( (ky_r-ky_real_r)./ky_real_r ).^2;
    end
    if abs(kxy_real_r)>1e-10
        Sum = Sum + ( (kxy_r-kxy_real_r)./kxy_real_r ).^2;
    end
    theta_error_all = sqrt(Sum);
    
    num_pred = find(theta_error_all==min(theta_error_all))
    num_pred = num_pred(1);
    theta_error = theta_error_all(num_pred);
end
%% error of a and b
a_error = abs((a_pred(num_pred)-a)/a);
b_error = abs((b_pred(num_pred)-b)/b);

end

