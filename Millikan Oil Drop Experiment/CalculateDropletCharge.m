%% CalculateDropletCharge.m
%
% Author: Daniel Opdahl
% Last modified: 11/10/2019
% Purpose: Takes user data from the Milikan Oil Drop experiment and
% processes it first by finding the radius and then the charge of a certain
% drop (along with uncertainties).

% Create arrays for charge values to be stored in

collected_charges = [];
collected_charges_unc = [];

% Manually inputted data
down_velos = [0.011010];
down_velos_unc = [0.000030];
up_velos = [0.036350];
up_velos_unc = [0.000147];
measured_viscosity = 1.8245;
measured_viscosity_unc = 0.01;
measured_voltage = 200;
measured_voltage_unc = 1;

% Use a for loop to run through all our manually inputted data

for i = 1:length(down_velos)
    
    % Input data, convert units, and define constants

    velocity_down = down_velos(i); %(mm/s)
    velocity_up = up_velos(i); %(mm/s)
    velocity_down_unc = down_velos_unc(i); %(mm/s)
    velocity_up_unc = up_velos_unc(i); %(mm/s)
    viscosity_air = measured_viscosity; %(Nsm^-2 * 10^-5)
    plate_separation = 0.00745; %(m)
    plate_separation_unc = 0.00001; %(m)
    voltage = measured_voltage; %(volts)

    velocity_down = velocity_down * 0.001; %(m/s)
    velocity_up = velocity_up * 0.001; %(m/s)
    velocity_down_unc = velocity_down_unc * 0.001; %(m/s)
    velocity_up_unc = velocity_up_unc * 0.001; %(m/s)
    viscosity_air = viscosity_air * 10^-5; %(Nsm^-2)
    measured_viscosity_unc = measured_viscosity_unc * 10^-5;

    density_oil = 866; %(kg/m^3)
    g = 9.81; %(m/s^2)

    % Calculate droplet radius

    droplet_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) );
    
    % Define partial derivatives for radius uncertainty
    
    dqDviscosity_air = (1/2) * ( (9*velocity_down*viscosity_air) / (2*density_oil*g) )^(-0.5) * ((9*velocity_down) / (2*density_oil*g));
    
    dqDvelocity_down = (1/2) * ( (9*velocity_down*viscosity_air) / (2*density_oil*g) )^(-0.5) * ((9*viscosity_air) / (2*density_oil*g));
    
    % Calculate uncertainty in droplet radius
    
    droplet_radius_unc = sqrt( (dqDviscosity_air * measured_viscosity_unc)^2 + (dqDvelocity_down * velocity_down_unc)^2 );

    % Calculate charge on droplet and add to collection

    charge = (6*pi*viscosity_air*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage);
    collected_charges(i) = charge;
    
    % Define partial derivatives for droplet charge uncertainty
    
    dqDviscosity_air = (6*pi*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage);
    
    dqDdroplet_radius = (6*pi*viscosity_air * (velocity_up + velocity_down) * plate_separation) / (voltage);
    
    dqDvelocity_up = (6*pi*viscosity_air*droplet_radius * plate_separation) / (voltage);
    
    dqDvelocity_down = (6*pi*viscosity_air*droplet_radius * plate_separation) / (voltage);
    
    dqDplate_separation = (6*pi*viscosity_air*droplet_radius * (velocity_up + velocity_down)) / (voltage);
    
    dqDvoltage = (-6*pi*viscosity_air*droplet_radius * (velocity_up + velocity_down) * plate_separation) / (voltage^2);
    
    % Calculate uncertainty in droplet charge and add to collection
    
    charge_unc = sqrt( (dqDviscosity_air*measured_viscosity_unc)^2 + (dqDdroplet_radius*droplet_radius_unc)^2 + (dqDvelocity_up*velocity_up_unc)^2 + (dqDvelocity_down*velocity_down_unc)^2 + (dqDplate_separation*plate_separation_unc)^2 + (dqDvoltage*measured_voltage_unc)^2 ) ;
    collected_charges_unc(i) = charge_unc;
end

% Display collections

collected_charges
collected_charges_unc