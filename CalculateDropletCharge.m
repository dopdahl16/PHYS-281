%% CalculateDropletCharge.m
%
% Author: Daniel Opdahl
% Last modified: 11/2/2019
% Purpose: Takes user input of a velocity of an oil drop and the viscosity of air
% in order to calculate the charge on the droplet for the Millikan Oil Drop
% Experiment.

% Calculate droplet radius
prompt = {'Input the measured velocity (mm/s)','Input the measured resistance (Megaohms)'};
dlg_title = 'Input parameters';
num_lines = 1;
answer = inputdlg(prompt, dlg_title, num_lines);
velocity = answer(1); %(mm/s)
resistance = answer(2); %(Megaohm)
velocity = velocity * 0.001; %(m/s)

density_oil = 886; %(kg/m^3)
g = 9.8; %(m/s^2)
sphere_radius = sqrt( (9*) / () )
