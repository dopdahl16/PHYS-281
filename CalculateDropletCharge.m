%% CalculateDropletCharge.m
%
% Author: Daniel Opdahl
% Last modified: 11/2/2019
% Purpose: Takes user input of Milikan Oil Drop experiment
% in order to calculate the charge on the droplet for the Millikan Oil Drop
% Experiment.

% Create arrays for charge values to be stored in

collected_charges = [];
collected_charges_unc = [];

% Manually inputted data
down_velos = [0.0187577];
down_velos_unc = [0.0000546];
up_velos = [0.0306450];
up_velos_unc = [0.0000954];
measured_viscosity = 1.819;
measured_viscosity_unc = 0.003;
measured_voltage = 200;
measured_voltage_unc = 1;

% Use a for loop to run through all our manually inputted data

for i = 1:length(down_velos)
    
    % Input data, convert units, and define constants

    velocity_down = down_velos(i); %(mm/s)
    velocity_up = up_velos(i); %(mm/s)
    viscosity_air = measured_viscosity; %(Nsm^-2 * 10^-5)
    plate_separation = 0.00745; %(m)
    plate_separation_unc = 0.00005; %(m)
    voltage = measured_voltage; %(volts)

    velocity_down = velocity_down * 0.001; %(m/s)
    velocity_up = velocity_up * 0.001; %(m/s)
    viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)
    measured_viscosity_unc = measured_viscosity_unc * 10^-5;

    density_oil = 886; %(kg/m^3)
    g = 9.81; %(m/s^2)

    % Calculate droplet radius

    droplet_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) )
    
    % Define partial derivatives for radius uncertainty
    
    dqDviscosity_air = (1/2) * sqrt( (9*velocity_down) / (2*density_oil*g) ) * (1/ sqrt(viscosity_air));
    
    dqDvelocity_down = (1/2) * sqrt( (9*viscosity_air) / (2*density_oil*g) ) * (1/ sqrt(velocity_down));
    
    % Calculate uncertainty in droplet radius
    
    droplet_radius_unc = sqrt( (dqDviscosity_air * measured_viscosity_unc)^2 + (dqDvelocity_down * down_velos_unc(i))^2 )

    % Calculate charge on droplet and add to collection

    charge = (6*pi*viscosity_air*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage)
    collected_charges(i) = charge;
    
    % Define partial derivatives for droplet charge uncertainty
    
    dqDviscosity_air = (6*pi*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage);
    
    dqDdroplet_radius = (6*pi*viscosity_air * (velocity_up + velocity_down) * plate_separation) / (voltage);
    
    dqDvelocity_up = (6*pi*viscosity_air*droplet_radius * plate_separation) / (voltage);
    
    dqDvelocity_down = (6*pi*viscosity_air*droplet_radius * plate_separation) / (voltage);
    
    dqDplate_separation = (6*pi*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage);
    
    dqDvoltage = (-6*pi*viscosity_air*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage)^2;
    
    % Calculate uncertainty in droplet charge
    
    charge_unc = sqrt( (dqDviscosity_air*measured_viscosity_unc)^2 + (dqDdroplet_radius*droplet_radius_unc)^2 + (dqDvelocity_up*down_velos_unc(i))^2 + (dqDvelocity_down*down_velos_unc(i))^2 + (dqDplate_separation*plate_separation_unc)^2 + (dqDvoltage*measured_voltage_unc)^2 )
    collected_charges_unc(i) = charge_unc;
end


    
% Take collected charges and display them in a histogram

% p=[1,1,1,2,4,5,5,5,5,5,8,10];
% h = histogram(p)
h = histogram(collected_charges)