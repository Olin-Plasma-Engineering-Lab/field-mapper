clear;
% field mapper measurement comparison script - Fall 2025 DR

% Background: In Spring 2025, we attempted to fire our HET, but ended up getting mixed results.
% We weren't able to get steady ionization (and as a result, a Hall Current). The primary focus of 
% the Fall 2025 semester has been to gather more data on why this might be the case. One facet of this
% analysis is a verification of the measurements of the field mapper against an expensive gaussmeter - 
% just to make sure that the field mapper is working as expected, and to eliminate it as a potential source
% of error.

% load in measured data
gaussmeter = readtable("gaussmeter.csv");
probe = readtable("mapper.csv");

% define plotting variables
gaussmeter_x = gaussmeter{:, 1}.*-1; 
gaussmeter_B = gaussmeter{:, 2};
probe_x = probe{:, 1}.*-1;
probe_B = probe{:, 2};

% plot data
figure;
plot(gaussmeter_x, gaussmeter_B, 'r--');
hold on;
grid on;
plot(gaussmeter_x, gaussmeter_B, 'ro');
xlabel("Axial Position from Exit Plane (mm)");
ylabel("Magnetic Field Strength (mT)");
title("Magnetic Field Strength (mT) vs. Axial Position (mm)");
subtitle("Comparison of Gaussmeter and Field Mapper Measurements");
plot(probe_x, probe_B, 'b--');
plot(probe_x, probe_B, 'bo');
legend('Gaussmeter', '', 'Field Mapper', '');
hold off;

%% statistical analysis

% compute percent error at each point
percent_error = abs((probe_B(1:4, :) - gaussmeter_B(1:4, :))./gaussmeter_B(1:4, :)) .* 100;
avg_percent_error = mean(percent_error);

sprintf("Average Percent Error: %.2f%%", avg_percent_error)

% plot results
figure;
plot(probe_x(1:4, :), percent_error, 'ro');
grid on;
xlabel("Axial Position from Exit Plane (mm)");
ylabel("Percent Error (%)");
title("Percent Error of Field Mapper vs. Gaussmeter");
hold off;