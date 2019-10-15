%% BubbleTrackDataAnalysis.m
%
% Author: Erin E. Flater, Luther College, Sept 2017
% Revised by: Daniel Opdahl
% Last modified: 9/24/19
%
% Purpose: Import, analyze, and plot electron bubble track data
%
%       Instructions: Follow the lab handout for details, but every line
%       which has an instruction in CAPITAL LETTERS should be replaced by a
%       variable or equation. Remember to finish every line of code with a
%       semicolon ; .
%%

%% Import data into MATLAB

TITLE= 'Select the Excel file that contains the data you want to bring into MATLAB';
[filename,filepath] = uigetfile('*.*', TITLE); %Prompts the user to select a data file
full_filename = fullfile( filepath, filename );

%Reading in the data from the file
data_matrix = xlsread(full_filename)

% Extract the first column of data as a vector
radius = data_matrix(:,2)

% Extract the second column of data as a vector
radius_unc = data_matrix(:,3)

% Extract the first column of data as a vector
distance = data_matrix(:,4)

% Extract the second column of data as a vector
distance_unc = data_matrix(:,5)


%% Rename variables
%The variables below should be the imported direct measurements, not yet converted to the scale of
% the paper.
r_meas = radius;                    %units: cm, data type: column vector
r_meas_unc = radius_unc;            %units: cm, data type: column vector
s_meas = distance;                  %units: cm, data type: column vector
s_meas_unc = distance_unc;          %units: cm, data type: column vector


%% Convert distances measured by a ruler/roller to units of the picture

conversion_factor = 0.707;

r_actual = conversion_factor * r_meas;
r_actual_unc = conversion_factor * r_meas_unc;
s_actual = conversion_factor * s_meas;
s_actual_unc = conversion_factor * s_meas_unc;

%% Use r and unc_r to calculate p (momentum) and unc_p (unc in momentum)
p = 2.08e-16*r_actual;                  %units: g*cm/s, data type: column vector
p_unc = 2.08e-16*r_actual_unc;          %units: g*cm/s, data type: column vector


%% Convert s and unc_s to KE (kinetic energy) and unc_KE (unc in KE)
slope = (0.04e-5);             %units: erg/cm, data type: scalar
intercept = 7e-7;              %units: erg, data type: scalar

KE = slope*s_actual + intercept;        %units: erg, data type: vector
KE_unc = slope*s_actual_unc;            %units: erg, data type: vector


%% Classical and relativistic KE
%mass of the electron in cgs units
m = 9.11e-31*1000;       %units: g
%speed of light in cgs units
c = 3e8*100;             %units: cm/s



%% Plot data with two theory curves: classical and relativitic theory

%Plot data with x and y error bars
errorbar(p,KE,KE_unc,KE_unc,p_unc,p_unc,'o')
hold on

%Plot theory curves with data
%Classical theory
p_class = 0:0.1e-16:5e-16;

KE_class = p_class.^2/(2*m);                  %units: erg, data type: vector
plot(p_class, KE_class)

%Relativistic theory
p_rel = 0:0.1e-16:20e-16;

KE_rel = sqrt((m*c^2)^2+(p_rel*c).^2)-(m*c^2);   %units: erg, data type: vector
plot(p_rel, KE_rel)

%ADD LABEL AXES
xlabel('momentum (cm*g/s)')
ylabel('KE (erg)')

hold off