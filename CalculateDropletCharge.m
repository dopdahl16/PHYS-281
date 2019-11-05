%% CalculateDropletCharge.m
%
% Author: Daniel Opdahl
% Last modified: 11/2/2019
% Purpose: Takes user input of a velocity of an oil drop and the viscosity of air
% in order to calculate the charge on the droplet for the Millikan Oil Drop
% Experiment.

% % % Calculate droplet radius

% % % prompt = {'Input the measured velocity down (mm/s)','Input the measured velocity up (mm/s)','Input the measured resistance converted to degrees Celcius'};
% % % dlg_title = 'Input parameters';
% % % num_lines = 1;
% % % answer = inputdlg(prompt, dlg_title, num_lines);
% % % velocity_down = answer(1); %(mm/s)
% % % velocity_up = answer(2); %(mm/s)
% % % temp = answer(3); %(C)


% Input data, convert units, and define constants

velocity_down = 0.187577; %(mm/s)
velocity_up = 0.306450; %(mm/s)
temp = 19; %(C)
plate_separation = 0.001; %(m)
voltage = 200; %(volts)

velocity_down = velocity_down * 0.001; %(m/s)
velocity_up = velocity_up * 0.001; %(m/s)

density_oil = 886; %(kg/m^3)
g = 9.81; %(m/s^2)

viscosity_air = 0.00475*temp + 1.7288; %(Nsm^-2 * 10^-5)
viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)

% Calculate droplet radius

sphere_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) )

% Calculate charge on droplet

charge = (6*pi*viscosity_air*sphere_radius* (velocity_up + velocity_down) *plate_separation) / (voltage)
