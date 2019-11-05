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


% Create array for charge values to be stored in

collected_charges = [];

% Manually inputted data for Tom
down_velos = [];
up_velos = [];
measured_temp = ;
measured_voltage = ;

% Use a for loop to run through all our manually inputted data

for i = 1:length(down_velos)
    
    % Input data, convert units, and define constants

    velocity_down = down_velos(i); %(mm/s)
    velocity_up = up_velos(i); %(mm/s)
    temp = measured_temp; %(C)
    plate_separation = 0.001; %(m)
    voltage = measured_voltage; %(volts)

    velocity_down = velocity_down * 0.001; %(m/s)
    velocity_up = velocity_up * 0.001; %(m/s)

    density_oil = 886; %(kg/m^3)
    g = 9.81; %(m/s^2)

    viscosity_air = 0.00475*temp + 1.7288; %(Nsm^-2 * 10^-5)
    viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)

    % Calculate droplet radius

    sphere_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) )

    % Calculate charge on droplet and add to collection

    charge = (6*pi*viscosity_air*sphere_radius* (velocity_up + velocity_down) *plate_separation) / (voltage)
    collected_charges.append(charge)
    
end

% Manually inputted data for Frankenstein
down_velos = [];
up_velos = [];
measured_temp = ;
measured_voltage = ;

% Use a for loop to run through all our manually inputted data

for i = 1:length(down_velos)
    
    % Input data, convert units, and define constants

    velocity_down = down_velos(i); %(mm/s)
    velocity_up = up_velos(i); %(mm/s)
    temp = measured_temp; %(C)
    plate_separation = 0.001; %(m)
    voltage = measured_voltage; %(volts)

    velocity_down = velocity_down * 0.001; %(m/s)
    velocity_up = velocity_up * 0.001; %(m/s)

    density_oil = 886; %(kg/m^3)
    g = 9.81; %(m/s^2)

    viscosity_air = 0.00475*temp + 1.7288; %(Nsm^-2 * 10^-5)
    viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)

    % Calculate droplet radius

    sphere_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) )

    % Calculate charge on droplet and add to collection

    charge = (6*pi*viscosity_air*sphere_radius* (velocity_up + velocity_down) *plate_separation) / (voltage)
    collected_charges.append(charge)
    
end

% Manually inputted data for Eric
down_velos = [];
up_velos = [];
measured_temp = ;
measured_voltage = ;

% Use a for loop to run through all our manually inputted data

for i = 1:length(down_velos)
    
    % Input data, convert units, and define constants

    velocity_down = down_velos(i); %(mm/s)
    velocity_up = up_velos(i); %(mm/s)
    temp = measured_temp; %(C)
    plate_separation = 0.001; %(m)
    voltage = measured_voltage; %(volts)

    velocity_down = velocity_down * 0.001; %(m/s)
    velocity_up = velocity_up * 0.001; %(m/s)

    density_oil = 886; %(kg/m^3)
    g = 9.81; %(m/s^2)

    viscosity_air = 0.00475*temp + 1.7288; %(Nsm^-2 * 10^-5)
    viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)

    % Calculate droplet radius

    sphere_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) )

    % Calculate charge on droplet and add to collection

    charge = (6*pi*viscosity_air*sphere_radius* (velocity_up + velocity_down) *plate_separation) / (voltage)
    collected_charges.append(charge)
    
end

% Manually inputted data for Mufasa
down_velos = [];
up_velos = [];
measured_temp = ;
measured_voltage = ;

% Use a for loop to run through all our manually inputted data

for i = 1:length(down_velos)
    
    % Input data, convert units, and define constants

    velocity_down = down_velos(i); %(mm/s)
    velocity_up = up_velos(i); %(mm/s)
    temp = measured_temp; %(C)
    plate_separation = 0.001; %(m)
    voltage = measured_voltage; %(volts)

    velocity_down = velocity_down * 0.001; %(m/s)
    velocity_up = velocity_up * 0.001; %(m/s)

    density_oil = 886; %(kg/m^3)
    g = 9.81; %(m/s^2)

    viscosity_air = 0.00475*temp + 1.7288; %(Nsm^-2 * 10^-5)
    viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)

    % Calculate droplet radius

    sphere_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) )

    % Calculate charge on droplet and add to collection

    charge = (6*pi*viscosity_air*sphere_radius* (velocity_up + velocity_down) *plate_separation) / (voltage)
    collected_charges.append(charge)
    
end
    
% Take collected charges and display them in a histogram

% p=[1,1,1,2,4,5,5,5,5,5,8,10];
% h = histogram(p)
h = histogram(collected_charges)