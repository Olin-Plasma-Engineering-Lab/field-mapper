% Validation file for the field mapper Hall Probe. This file has two
% different possible plots: one for the 130mT chip, and one for a combined
% 5mT chip (at lower field strengths) and 130mT chip (at higher field
% strengths). The file also notes the MAPE to try and quantify the error
% between the probe measurements and reality.
clear;

% load probe data
probe = readtable(""); % ENTER FILE NAME

% % Compute predicted B-field strengths (T)
% % UNCOMMENT FOR 130mT CHIP
% 
% u_0 = 4*pi*(10e-7); % in T-m/A
% m = 0.83; % in T - magnet strength
% l = 0.0047625; % in m
% positions = [2; 5; 7; 10; 12; 15; 17; 20]./1000 + l; % probe positions, in m
% M = m*2*l; % magnetic moment
% predicted_b = ((u_0/(4*pi)) * ((2*M)./(positions.^3))) * 1000; % in mT
% 
% % Compute normalized B-field strengths (mT)
% 
% measured_b = [];
% for i = 1:length(probe.DataSet1_XMagneticField130mT_mT_)
%     current_b = [probe.DataSet1_XMagneticField130mT_mT_(i); 
%         probe.DataSet1_YMagneticField130mT_mT_(i); 
%         probe.DataSet1_ZMagneticField130mT_mT_(i)];
%     measured_b = [measured_b; norm(current_b)];
% end


% Compute predicted B-field strengths (T)
% UNCOMMENT FOR SENSOR COMBINATION

u_0 = 4*pi*(10e-7); % in T-m/A
m = 0.83; % in T - magnet strength
l = 0.0047625; % in m
positions = [2; 5; 7; 10; 12; 15; 17; 20]./1000 + l; % probe positions, in m
M = m*2*l; % magnetic moment
predicted_b = ((u_0/(4*pi)) * ((2*M)./(positions.^3))) * 1000; % in mT

% Compute normalized B-field strengths (mT)
measured_b = [];
for i = 1:length(probe.DataSet1_XMagneticField130mT_mT_)
    % Use 5mT chip measurements below 10ish mT
    if i > 2
        current_b = [probe.DataSet1_XMagneticField_mT_(i); 
            probe.DataSet1_YMagneticField_mT_(i); 
            probe.DataSet1_ZMagneticField_mT_(i);];
        measured_b = [measured_b; norm(current_b)];
    else 
        current_b = [probe.DataSet1_XMagneticField130mT_mT_(i); 
            probe.DataSet1_YMagneticField130mT_mT_(i); 
            probe.DataSet1_ZMagneticField130mT_mT_(i)];
        measured_b = [measured_b; norm(current_b)];
    end
end

% Compare predicted and measured curves
figure;
plot(positions, measured_b, 'r', 'LineWidth', 2)
hold on;
grid on;
plot(positions, predicted_b, 'b', 'LineWidth', 2)
title('Hall Probe Error Characterization')
subtitle('Predicted vs. Measured B-field strength')
xlabel('Distance (m)')
ylabel('Field Strength (mT)')
legend('Measured Strength', 'Predicted Strength')


% Error Quantification
error = mape(predicted_b, measured_b);
display = sprintf('MAPE is %.4f%%', error);
disp(display)