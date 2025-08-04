% Data visualization script for the magnetic field mapper (Scitech 2025)
clear;

% Continuous Data Template
% The structure of this file is as follows: lines 1-67 do all of the
% computation, including processing the raw csv data, combined mapping and
% positions data, and ultimately downsampling the resulting matrix to
% increase plotting resolution. The rest of the script plots the data in
% different ways to get a more holistic picture of what is going on inside
% the HET channel. Each plot has its own figure window.

% Make sure that the mapping data csv and COMSOL 1D simulation .txt file
% are both in the same directory as this script. Any place with a comment
% in all caps indicates that a parameter should be added/edited.

% load continuous data
mapping_data = struct2cell(data_process());

% mapping parameters
column_number = ; % CHANGE BASED ON NUMBER OF MAPPING COLUMNS
col_dist = ; % in mm; CHANGE BASED ON THE DISTANCE BETWEEN MAPPING COLUMNS
speed = 50/60; % extruder speed in mm/s, F50; CHANGE BASED ON EXTRUDER SPEED
first_point = [0; 0; 0];
wall_offset = ; % in mm; CHANGE BASED ON PROBE OFFSET FROM THE WALL (DON'T ACCOUNT FOR PROBE THICKNESS)
channel_width = ; % in mm; CHANGE BASED ON THRUSTER DIMENSIONS
channel_depth = ; % in mm; CHANGE BASED ON THRUSTER DIMENSIONS
% distance from chip - chip length/thickness (all from datasheet) -> measurements taken in center
% of chip
dist_from_chipx = 1.76 - 0.55;
dist_from_chipz = 10.5 - 1.45;

% x-coordinates - no offsets
x_coordinates = zeros(1, column_number);
next_point = 0;
for i = 1:column_number
    x_coordinates(i) = first_point(1) + next_point;
    next_point = next_point + col_dist;
end

% create new cells containing x-pos, z-pos, and field data
b_field_data = {};
positions_data = {};
for i = 1:column_number
    % positions - account for negative movement of the probe into the channel + probe
    % offset (10.5 - 1.45 mm in z, 1.76 - 0.55 mm in x); center of chip for z,
    % x (thickness of sensor is 1.1 mm)

    % Make position data dimensionless by dividing by channel width and
    % depth (x and z respectively)    
    z_positions = ((mapping_data{i, 1}(:, 1) * speed * -1) + dist_from_chipz)./channel_depth;
    x_positions = (zeros(length(z_positions), 1) + x_coordinates(i) + dist_from_chipx + wall_offset)./channel_width;
    
    % final data placement
    positions_data{i} = [x_positions, z_positions];
end

for i = 1:column_number
    % b-field data - normalized vectors, and magnitude calcs (4 columns)
    temp_b = zeros(length(mapping_data{i, 1}(:, 5:7)), 4);
    for j = 1:length(mapping_data{i, 1}(:, 5:7))
        temp_b(j, :) = b_field_process(mapping_data{i, 1}(j, 5:7));
    end
    b_field_data{i} = temp_b;

end

% create massive matrix for plotting (x, z, bx, by, bz, bmag)
data_matrix = [];
for i = 1:column_number
    temp = [positions_data{1, i}, b_field_data{1, i}];
    data_matrix = [data_matrix; temp];
end

% Cut down dataset to increase plotting resolution (downsampling)
data_matrix_cleaned = data_matrix(1:5:14900,:);

% COMSOL 1D data to overlay with measured results
comsol_1D = readmatrix(""); % EDIT WITH FILE NAME
% transform raw COMSOL data to the ideal reference frame (exit plane = 0)
mag_sim_z = (comsol_1D(:, 1) * -1 * 1000 + 10)/channel_depth; % EDIT IF NEEDED, BASED ON RAW DATA

% curve fitting to get a continuous polynomial to represent the simulation
fitted_z = polyfit(mag_sim_z, comsol_1D(:, 2), 7);
fitted_model = polyval(fitted_z, mag_sim_z);
% UNCOMMENT LINE BELOW TO TEST WHICH POLYNOMIAL ORDER FITS BEST
%residual = mape(comsol_1D(:, 2), fitted_model)

%PLOTTING%

% 1D B-field plot (magnitude); center of channel
figure;
% Compute the closest column number at the thruster centerline; based on
% screen dimensions
mid_column_position = ; % in mm; CHANGE BASED ON THRUSTER DIMENSIONS
mid_column = round((mid_column_position + (dist_from_chipx * -1) + wall_offset)/col_dist);
% plot simulated data, fitted curve, and measured data
plot(positions_data{1, mid_column}(:, 2), b_field_data{1, mid_column}(:, 4), 'bo', 'MarkerSize', 1)
hold on;
plot(mag_sim_z, comsol_1D(:, 2), 'ro', 'MarkerSize', 4)
plot(mag_sim_z, fitted_model, 'k')
title('1D Magnitude Plot')
subtitle('B-Field Strength vs. Z-position from the Exit Plane')
xlabel('Position from the Exit Plane (Channel Depth Ratio)')
ylabel('Magnetic Field Strength (mT)')
xlim([-1 1])
ylim([0 30])
legend('Measured', 'Simulated', 'Fitted')
hold off;

% Compute residuals for 1D magnitude plot
figure;
z_pos_channel_depth = positions_data{1, mid_column}(:, 2);
fit_to_measurement_mT = polyval(fitted_z, z_pos_channel_depth);
residuals = b_field_data{1, mid_column}(:, 4) - fit_to_measurement_mT;
plot(z_pos_channel_depth, residuals, 'o');
title('1D Magnitude Residuals Plot')
subtitle('Residuals vs. Z-position from the Exit Plane')
xlabel('Position from the Exit Plane (Channel Depth Ratio)')
ylabel('B-field Strength Residuals (mT)')
hold off;

% 2D vector plot - shows direction of the B-field at points in the HET
% channel
figure;
quiver(data_matrix_cleaned(:, 1), data_matrix_cleaned(:, 2), data_matrix_cleaned(:, 3), data_matrix_cleaned(:, 5), 'LineWidth', 0.5)
title('2D Directionality Plot - HET Channel Cross Section')
subtitle('Z-position vs. X-position')
grid on;
axis equal;
xlabel('X-position (Channel Width Ratio)')
ylabel('Z-position (Channel Depth Ratio)')
hold off;

% 2D Streamline plot - shows magnetic field streamlines
figure;
% Meshgrid of unique X and Z values
[X, Z] = meshgrid(unique(data_matrix_cleaned(:, 1)), unique(data_matrix_cleaned(:, 2)));
% Interpolation
U = griddata(data_matrix_cleaned(:, 1), data_matrix_cleaned(:, 2), data_matrix_cleaned(:, 3), X, Z);
V = griddata(data_matrix_cleaned(:, 1), data_matrix_cleaned(:, 2), data_matrix_cleaned(:, 5), X, Z);

% Starting points - meshgrid; nxn grid
grid_resolution = 12; % CHANGE TO INCREASE PLOT RESOLUTION
[startX, startY] = meshgrid(linspace(min(data_matrix_cleaned(:, 1)), max(data_matrix_cleaned(:, 1)), grid_resolution), ...
    linspace(min(data_matrix_cleaned(:, 2)), max(data_matrix_cleaned(:, 2)), grid_resolution));
% Starting points - along Y axis
% startY = linspace(min(data_matrix_cleaned(:, 2)), max(data_matrix_cleaned(:, 2)), 100);
% startX = zeros(size(startY)) + min(data_matrix_cleaned(:, 1);
streamline(X, Z, U, V, startX, startY)
title('2D Streamline Plot - HET Channel Cross Section')
subtitle('Z-position vs. X-position')
grid on;
axis equal;
xlabel('X-position (Channel Length Ratio)')
ylabel('Z-position (Channel Depth Ratio)')
hold off;

% 2D Strength Contours plot - shows magnetic field strength contours
figure;
% Interpolation
field_strength = griddata(data_matrix(:, 1), data_matrix(:, 2), data_matrix(:, 6), X, Z, 'cubic');
contourf(X, Z, field_strength);
c = colorbar;
c.Label.String = 'B-field Strength (mT)';
colormap hot;
title('2D Strength Contour Plot - HET Channel Cross Section')
subtitle('Z-position vs. X-position')
grid on;
axis equal;
xlabel('X-position (Channel Width Ratio)')
ylabel('Z-position (Channel Depth Ratio)')
hold off;

% Mock Comsol Plot - Strength + Streamlines, w/channel and screen walls
figure;
% Plot strength contours in the background
contourf(X, Z, field_strength)
c = colorbar;
c.Label.String = 'B-field Strength (mT)';
colormap hot;
% layer streamlines
h = streamline(X, Z, U, V, startX, startY);
hold on;
set(h, 'Color', [0 1 1])
title('2D Strength Contour + Streamlines Plot')
subtitle('Z-position vs. X-position, 0 degrees')
grid on;
%axis equal;
xlabel('X-position (Channel Width Ratio)')
ylabel('Z-position (Channel Depth Ratio)')

% Plot Channel and Screen Walls as lines (Dimensions in Solidworks)
x_screen = [0 0 15.12 15.12 17.12 17.12 -2.04 -2.04 0]./channel_width;
y_screen = [-3.5 -18.5 -18.5 -3.5 -3.5 -20.5 -20.5 -3.5 -3.5]./channel_depth;
x_channel = [0 2 2 13.08 13.08 15.12 15.12 0 0]./channel_width;
y_channel = [0 -3.46 -16.5 -16.5 -3.46 0 -18.5 -18.5 0]./channel_depth;
% add relevant geometry
plot(x_screen, y_screen, 'k', 'LineWidth', 1.5)
plot(x_channel, y_channel, 'r', 'LineWidth', 1.5)
set(gca, 'color', [0.88 0.88 0.88]);
legend_cust = [plot(nan, 'c'), plot(nan, 'r'), plot(nan, 'k')];
legend(legend_cust, {'Streamlines', 'HET Channel', 'HET Screen'}, 'location', 'northeastoutside');
hold off;


% DATA PROCESSING FUNCTIONS

function data = data_process()
    data_matrix = readmatrix(""); % RENAME BASED ON MAPPING DATA FILE
    data = [];
    % generate a structure to split complete dataset into individual column
    % datasets
    for i = 1:size(data_matrix, 2)/7
        temp_cell = data_matrix(:, ((i-1)*6+i):(i*6+i));
        temp_cell = rmmissing(temp_cell);
        data.("col"+i) = temp_cell;
    end
end

% Take raw B-field data, find the magnitude, and normalize the vector
function b_field_processed = b_field_process(field_vector)
    % convert input to reference frame of the probe
    field_vector = [field_vector(3), field_vector(2), field_vector(1)*-1];
    b_mag = norm(field_vector);
    b_norm = field_vector./b_mag;
    b_field_processed = [b_norm(1), b_norm(2), b_norm(3), b_mag];
end


