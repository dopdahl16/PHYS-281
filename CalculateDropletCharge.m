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
down_velos = [0.0146,
0.02032,
0.020585,
0.017322,
0.018723,
0.016018,
0.017739];
down_velos_unc = [0.00004,
0.0001,
0.00007,
0.00005,
0.00004,
0.00004,
0.00007];
up_velos = [0.1524,
0.6803,
0.8675,
0.4535,
0.5337,
0.5801,
0.15261];
up_velos_unc = [0.0004,
0.003,
0.003,
0.001,
0.001,
0.002,
0.0006];
measured_viscosity = 1.819;
measured_viscosity_unc = 0.1;
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

    density_oil = 866; %(kg/m^3)
    g = 9.81; %(m/s^2)

    % Calculate droplet radius

    droplet_radius = sqrt( (9*viscosity_air*velocity_down) / (2*density_oil*g) );
    
    % Define partial derivatives for radius uncertainty
    
    dqDviscosity_air = (1/2) * ( (9*velocity_down*viscosity_air) / (2*density_oil*g) )^(-0.5) * ((9*velocity_down) / (2*density_oil*g));
    
    dqDvelocity_down = (1/2) * ( (9*velocity_down*viscosity_air) / (2*density_oil*g) )^(-0.5) * ((9*viscosity_air) / (2*density_oil*g));
    
    % Calculate uncertainty in droplet radius
    
    droplet_radius_unc = sqrt( (dqDviscosity_air * measured_viscosity_unc)^2 + (dqDvelocity_down * (down_velos_unc(i)/1000))^2 );

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
    
    % Calculate uncertainty in droplet charge
    
    charge_unc = sqrt( (dqDviscosity_air*measured_viscosity_unc)^2 + (dqDdroplet_radius*droplet_radius_unc)^2 + (dqDvelocity_up*(up_velos_unc(i)/1000))^2 + (dqDvelocity_down*(down_velos_unc(i)/1000))^2 + (dqDplate_separation*plate_separation_unc)^2 + (dqDvoltage*measured_voltage_unc)^2 );
    collected_charges_unc(i) = charge_unc;
end

collected_charges
collected_charges_unc
    
% Take collected charges and display them in a histogram

% p=[1,1,1,2,4,5,5,5,5,5,8,10];
% h = histogram(p)
% h = histogram(collected_charges, 'BinEdges', linspace(0, 1, 5))
% h = histogram(collected_charges,10000000)
% grid on
% xlim([0,0.00000000000000000001]);

drop_number = 1:7;
charge = [8.00,39.59,50.51,24.57,29.97,29.91,8.99];
unc_charge = [0.4967,0.3848,0.4639,0.2151,0.2586,0.2723,0.0837];



x1 = drop_number;
y1 = charge;
%dx1 = errangle;
dy1 = unc_charge;

%Plot data with x and y error bars
errorbar(x1,y1,dy1,'bo')

grid on

xlabel('oil drop number')
ylabel('Charge x 10^-19 (C)')

% spacing = .1602;
% yticks(spacing*(1:100))

axis([1 7 0 50])

hold off