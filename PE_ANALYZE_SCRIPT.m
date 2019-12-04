
% Import Data
TITLE= 'Select the file with the data you want to bring into MATLAB';
[filename,filepath] = uigetfile('*.*', TITLE); %Prompts the user to select a data file
full_filename = fullfile( filepath, filename );

[~,SheetNames] = xlsfinfo(full_filename);
nSheets = length(SheetNames);
Data = [];
for ii=1:nSheets
    Name =  SheetNames{ii};
    Data = [Data, xlsread(full_filename, Name)];
end


forward_voltage_V577 = Data(:,1);
forward_current_A577 = Data(:,2);
reverse_voltage_V577 = Data(:,4);
reverse_current_A577 = Data(:,5);

forward_voltage_V546 = Data(:,6);
forward_current_A546 = Data(:,7);
reverse_voltage_V546 = Data(:,9);
reverse_current_A546 = Data(:,10);

forward_voltage_V436 = Data(:,11);
forward_current_A436 = Data(:,12);
reverse_voltage_V436 = Data(:,14);
reverse_current_A436 = Data(:,15);

forward_voltage_V405 = Data(:,16);
forward_current_A405 = Data(:,17);
reverse_voltage_V405 = Data(:,19);
reverse_current_A405 = Data(:,20);

forward_voltage_V365 = Data(:,21);
forward_current_A365 = Data(:,22);
reverse_voltage_V365 = Data(:,24);
reverse_current_A365 = Data(:,25);









% % 577 analysis
% 
% %%% FORWARD
% % Calculate average of three runs
% firstrun577 = forward_current_A577(1:301);
% secondrun577 = forward_current_A577(302:602);
% thirdrun577 = forward_current_A577(603:903);
% 
% averagecurrent577 = [];
% current_standard_dev_577 = [];
% for i = 1:301;
%     averagecurrent577i = (firstrun577(i) + secondrun577(i) + thirdrun577(i))/3;
%     current_standard_dev_577(i) = std([firstrun577(i); secondrun577(i); thirdrun577(i)]);
%     averagecurrent577(i) =  averagecurrent577i;
% end
% voltagetouse577 = forward_voltage_V577(1:301);
% 
% % Uncertainty in Voltage
% voltage_unc_577 = 0.001*ones(size(voltagetouse577));
% 
% % Uncertainty in Current
% current_standard_error_577 = current_standard_dev_577/ sqrt(3);
% 
% 
% 
% 
% 
% %%% REVERSE
% % Calculate average of three runs 
% firstrun_reverse_577 = reverse_current_A577(1:41);
% secondrun_reverse_577 = reverse_current_A577(42:82);
% thirdrun_reverse_577 = reverse_current_A577(83:123);
% 
% averagecurrent_reverse_577 =[];
% current_reverse_standard_dev_577 = [];
% for i = 1:41;
%     averagecurrent_reverse_577i = (firstrun_reverse_577(i) + secondrun_reverse_577(i) + thirdrun_reverse_577(i)) / 3;
%     current_reverse_standard_dev_577(i) = std([firstrun577(i); secondrun577(i); thirdrun577(i)]);
%     averagecurrent_reverse_577(i) = averagecurrent_reverse_577i;
% end
% voltagetouse_reverse_577 = reverse_voltage_V577(1:41);
% voltagetouse_reverse_577 = -voltagetouse_reverse_577;
% 
% % Uncertainty in Voltage
% voltage_unc_reverse_577 = 0.001*ones(size(voltagetouse_reverse_577));
% 
% % Uncertainty in Current
% current_reverse_standard_error_577 = current_reverse_standard_dev_577/ sqrt(3);
% 
% 
% 
% 
% % Plot a figure with both data sets on one graph
% figure
% errorbar(voltagetouse577, averagecurrent577',current_standard_error_577,current_standard_error_577,voltage_unc_577,voltage_unc_577);
% hold on;
% errorbar(voltagetouse_reverse_577, averagecurrent_reverse_577,current_reverse_standard_error_577,current_reverse_standard_error_577,voltage_unc_reverse_577,voltage_unc_reverse_577)
% plot(voltagetouse577, averagecurrent577');
% plot(voltagetouse_reverse_577, averagecurrent_reverse_577);
% title("Current vs. Voltage (577nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% 
% 
% % Plot just the forward bias
% figure
% errorbar(voltagetouse577, averagecurrent577',current_standard_error_577,current_standard_error_577,voltage_unc_577,voltage_unc_577);
% title("Forward Bias Current vs. Voltage (577nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% % Plot just the reverse bias 
% figure
% errorbar(voltagetouse_reverse_577, averagecurrent_reverse_577,current_reverse_standard_error_577,current_reverse_standard_error_577,voltage_unc_reverse_577,voltage_unc_reverse_577)
% title("Reverse Bias Current vs. Voltage (577nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% 
% 
% % Find knee using plateau uncertaities and deviation from that value
% 
% % Calculate aveage data point
% running_total_577 = 0;
% g577 = current_reverse_standard_error_577(35:40);
% for i = 1:6;
%     running_total_577 = g577(i) + running_total_577;
% end
% average_threshold_point_577 = running_total_577 / 6;
% 
% % Calculate average error bar size
% running_total_577 = 0;
% for i = 1:6;
%     running_total_577 = g577(i)^2 + running_total_577;
% end
% uncertainty_mean_577 = sqrt(running_total_577)/6;
% 
% % Calculate standard deviation
% standard_dev_577 = std(current_reverse_standard_error_577(35:40));
% 
% % Calculate threshold
% threshold_bar_577 = sqrt(uncertainty_mean_577^2 + standard_dev_577^2);
% 
% % Go through each data point after index 35, going down and if it plus its
% % error bar is out of the threshold range, then that point is the knee% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 546 analysis
% 
% %%% FORWARD
% % Calculate average of three runs
% firstrun546 = forward_current_A546(1:301);
% secondrun546 = forward_current_A546(302:602);
% thirdrun546 = forward_current_A546(603:903);
% 
% averagecurrent546 = [];
% current_standard_dev_546 = [];
% for i = 1:301;
%     averagecurrent546i = (firstrun546(i) + secondrun546(i) + thirdrun546(i))/3;
%     current_standard_dev_546(i) = std([firstrun546(i); secondrun546(i); thirdrun546(i)]);
%     averagecurrent546(i) =  averagecurrent546i;
% end
% voltagetouse546 = forward_voltage_V546(1:301);
% 
% % Uncertainty in Voltage
% voltage_unc_546 = 0.001*ones(size(voltagetouse546));
% 
% % Uncertainty in Current
% current_standard_error_546 = current_standard_dev_546/ sqrt(3);
% 
% 
% 
% 
% 
% %%% REVERSE
% % Calculate average of three runs 
% firstrun_reverse_546 = reverse_current_A546(1:41);
% secondrun_reverse_546 = reverse_current_A546(42:82);
% thirdrun_reverse_546 = reverse_current_A546(83:123);
% 
% averagecurrent_reverse_546 =[];
% current_reverse_standard_dev_546 = [];
% for i = 1:41;
%     averagecurrent_reverse_546i = (firstrun_reverse_546(i) + secondrun_reverse_546(i) + thirdrun_reverse_546(i)) / 3;
%     current_reverse_standard_dev_546(i) = std([firstrun546(i); secondrun546(i); thirdrun546(i)]);
%     averagecurrent_reverse_546(i) = averagecurrent_reverse_546i;
% end
% voltagetouse_reverse_546 = reverse_voltage_V546(1:41);
% voltagetouse_reverse_546 = -voltagetouse_reverse_546;
% 
% % Uncertainty in Voltage
% voltage_unc_reverse_546 = 0.001*ones(size(voltagetouse_reverse_546));
% 
% % Uncertainty in Current
% current_reverse_standard_error_546 = current_reverse_standard_dev_546/ sqrt(3);
% 
% 
% 
% 
% % Plot a figure with both data sets on one graph
% figure
% errorbar(voltagetouse546, averagecurrent546',current_standard_error_546,current_standard_error_546,voltage_unc_546,voltage_unc_546);
% hold on;
% errorbar(voltagetouse_reverse_546, averagecurrent_reverse_546,current_reverse_standard_error_546,current_reverse_standard_error_546,voltage_unc_reverse_546,voltage_unc_reverse_546)
% plot(voltagetouse546, averagecurrent546');
% plot(voltagetouse_reverse_546, averagecurrent_reverse_546);
% title("Current vs. Voltage (546nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% 
% 
% % Plot just the forward bias
% figure
% errorbar(voltagetouse546, averagecurrent546',current_standard_error_546,current_standard_error_546,voltage_unc_546,voltage_unc_546);
% title("Forward Bias Current vs. Voltage (546nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% % Plot just the reverse bias 
% figure
% errorbar(voltagetouse_reverse_546, averagecurrent_reverse_546,current_reverse_standard_error_546,current_reverse_standard_error_546,voltage_unc_reverse_546,voltage_unc_reverse_546)
% title("Reverse Bias Current vs. Voltage (546nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% 
% 
% % Find knee using plateau uncertaities and deviation from that value
% 
% % Calculate aveage data point
% running_total_546 = 0;
% g546 = current_reverse_standard_error_546(35:40);
% for i = 1:6;
%     running_total_546 = g546(i) + running_total_546;
% end
% average_threshold_point_546 = running_total_546 / 6;
% 
% % Calculate average error bar size
% running_total_546 = 0;
% for i = 1:6;
%     running_total_546 = g546(i)^2 + running_total_546;
% end
% uncertainty_mean_546 = sqrt(running_total_546)/6;
% 
% % Calculate standard deviation
% standard_dev_546 = std(current_reverse_standard_error_546(35:40));
% 
% % Calculate threshold
% threshold_546 = sqrt(uncertainty_mean_546^2 + standard_dev_546^2)
% 
% % Go through each data point after index 35, going down and if it plus its
% % error bar is out of the threshold range, then that point is the knee
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 436 analysis
% 
% %%% FORWARD
% % Calculate average of three runs
% firstrun436 = forward_current_A436(1:301);
% secondrun436 = forward_current_A436(302:602);
% thirdrun436 = forward_current_A436(603:903);
% 
% averagecurrent436 = [];
% current_standard_dev_436 = [];
% for i = 1:301;
%     averagecurrent436i = (firstrun436(i) + secondrun436(i) + thirdrun436(i))/3;
%     current_standard_dev_436(i) = std([firstrun436(i); secondrun436(i); thirdrun436(i)]);
%     averagecurrent436(i) =  averagecurrent436i;
% end
% voltagetouse436 = forward_voltage_V436(1:301);
% 
% % Uncertainty in Voltage
% voltage_unc_436 = 0.001*ones(size(voltagetouse436));
% 
% % Uncertainty in Current
% current_standard_error_436 = current_standard_dev_436/ sqrt(3);
% 
% 
% 
% %%% REVERSE
% % Calculate average of three runs 
% firstrun_reverse_436 = reverse_current_A436(1:41);
% secondrun_reverse_436 = reverse_current_A436(42:82);
% thirdrun_reverse_436 = reverse_current_A436(83:123);
% 
% averagecurrent_reverse_436 =[];
% current_reverse_standard_dev_436 = [];
% for i = 1:41;
%     averagecurrent_reverse_436i = (firstrun_reverse_436(i) + secondrun_reverse_436(i) + thirdrun_reverse_436(i)) / 3;
%     current_reverse_standard_dev_436(i) = std([firstrun436(i); secondrun436(i); thirdrun436(i)]);
%     averagecurrent_reverse_436(i) = averagecurrent_reverse_436i;
% end
% voltagetouse_reverse_436 = reverse_voltage_V436(1:41);
% voltagetouse_reverse_436 = -voltagetouse_reverse_436;
% 
% % Uncertainty in Voltage
% voltage_unc_reverse_436 = 0.001*ones(size(voltagetouse_reverse_436));
% 
% % Uncertainty in Current
% current_reverse_standard_error_436 = current_reverse_standard_dev_436/ sqrt(3);
% 
% 
% 
% 
% % Plot a figure with both data sets on one graph
% figure
% errorbar(voltagetouse436, averagecurrent436',current_standard_error_436,current_standard_error_436,voltage_unc_436,voltage_unc_436);
% hold on;
% errorbar(voltagetouse_reverse_436, averagecurrent_reverse_436,current_reverse_standard_error_436,current_reverse_standard_error_436,voltage_unc_reverse_436,voltage_unc_reverse_436)
% plot(voltagetouse436, averagecurrent436');
% plot(voltagetouse_reverse_436, averagecurrent_reverse_436);
% title("Current vs. Voltage (436nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% 
% 
% % Plot just the forward bias
% figure
% errorbar(voltagetouse436, averagecurrent436',current_standard_error_436,current_standard_error_436,voltage_unc_436,voltage_unc_436);
% title("Forward Bias Current vs. Voltage (436nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% % Plot just the reverse bias 
% figure
% errorbar(voltagetouse_reverse_436, averagecurrent_reverse_436,current_reverse_standard_error_436,current_reverse_standard_error_436,voltage_unc_reverse_436,voltage_unc_reverse_436)
% title("Reverse Bias Current vs. Voltage (436nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% 
% 
% % Find knee using plateau uncertaities and deviation from that value
% 
% % Calculate aveage data point
% running_total_436 = 0;
% g436 = current_reverse_standard_error_436(35:40);
% for i = 1:6;
%     running_total_436 = g436(i) + running_total_436;
% end
% average_threshold_point_436 = running_total_436 / 6;
% 
% % Calculate average error bar size
% running_total_436 = 0;
% for i = 1:6;
%     running_total_436 = g436(i)^2 + running_total_436;
% end
% uncertainty_mean_436 = sqrt(running_total_436)/6;
% 
% % Calculate standard deviation
% standard_dev_436 = std(current_reverse_standard_error_436(35:40));
% 
% % Calculate threshold
% threshold_436 = sqrt(uncertainty_mean_436^2 + standard_dev_436^2)
% 
% % Go through each data point after index 35, going down and if it plus its
% % error bar is out of the threshold range, then that point is the knee
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 405 analysis
% 
% %%% FORWARD
% % Calculate average of three runs
% firstrun405 = forward_current_A405(1:301);
% secondrun405 = forward_current_A405(302:602);
% thirdrun405 = forward_current_A405(603:903);
% 
% averagecurrent405 = [];
% current_standard_dev_405 = [];
% for i = 1:301;
%     averagecurrent365i = (firstrun405(i) + secondrun405(i) + thirdrun405(i))/3;
%     current_standard_dev_405(i) = std([firstrun405(i); secondrun405(i); thirdrun405(i)]);
%     averagecurrent405(i) =  averagecurrent405i;
% end
% voltagetouse405 = forward_voltage_V405(1:301);
% 
% % Uncertainty in Voltage
% voltage_unc_405 = 0.001*ones(size(voltagetouse405));
% 
% % Uncertainty in Current
% current_standard_error_405 = current_standard_dev_405/ sqrt(3);
% 
% 
% 
% 
% %%%REVERSE
% % Calculate average of three runs 
% firstrun_reverse_405 = reverse_current_A405(1:41);
% secondrun_reverse_405 = reverse_current_A405(42:82);
% thirdrun_reverse_405 = reverse_current_A405(83:123);
% 
% averagecurrent_reverse_405 =[];
% current_reverse_standard_dev_405 = [];
% for i = 1:41;
%     averagecurrent_reverse_405i = (firstrun_reverse_405(i) + secondrun_reverse_405(i) + thirdrun_reverse_405(i)) / 3;
%     current_reverse_standard_dev_405(i) = std([firstrun405(i); secondrun405(i); thirdrun405(i)]);
%     averagecurrent_reverse_405(i) = averagecurrent_reverse_405i;
% end
% voltagetouse_reverse_405 = reverse_voltage_V405(1:41);
% voltagetouse_reverse_405 = -voltagetouse_reverse_405;
% 
% % Uncertainty in Voltage
% voltage_unc_reverse_405 = 0.001*ones(size(voltagetouse_reverse_405));
% 
% % Uncertainty in Current
% current_reverse_standard_error_405 = current_reverse_standard_dev_405/ sqrt(3);
% 
% 
% % Plot a figure with both data sets on one graph
% figure
% errorbar(voltagetouse405, averagecurrent405',current_standard_error_405,current_standard_error_405,voltage_unc_405,voltage_unc_405);
% hold on;
% errorbar(voltagetouse_reverse_405, averagecurrent_reverse_405,current_reverse_standard_error_405,current_reverse_standard_error_405,voltage_unc_reverse_405,voltage_unc_reverse_405)
% plot(voltagetouse405, averagecurrent405');
% plot(voltagetouse_reverse_405, averagecurrent_reverse_405);
% title("Current vs. Voltage (365nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% 
% % Plot just the forward bias
% figure
% errorbar(voltagetouse405, averagecurrent405',current_standard_error_405,current_standard_error_405,voltage_unc_405,voltage_unc_405);
% title("Forward Bias Current vs. Voltage (405nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% % Plot just the reverse bias 
% figure
% errorbar(voltagetouse_reverse_405, averagecurrent_reverse_405,current_reverse_standard_error_405,current_reverse_standard_error_405,voltage_unc_reverse_405,voltage_unc_reverse_405)
% title("Reverse Bias Current vs. Voltage (405nm Filter)")
% xlabel('Voltage (V)') 
% ylabel('Current (A)')
% hold on;
% 
% 
% % Find knee using plateau uncertaities and deviation from that value
% 
% % Calculate aveage data point
% running_total_405 = 0;
% g405 = current_reverse_standard_error_405(35:40);
% for i = 1:6;
%     running_total_405 = g405(i) + running_total_405;
% end
% average_threshold_point_405 = running_total_405 / 6;
% 
% % Calculate average error bar size
% running_total_405 = 0;
% for i = 1:6;
%     running_total_405 = g405(i)^2 + running_total_405;
% end
% uncertainty_mean_405 = sqrt(running_total_405)/6;
% 
% % Calculate standard deviation
% standard_dev_405 = std(current_reverse_standard_error_405(35:40));
% 
% % Calculate threshold
% threshold_405 = sqrt(uncertainty_mean_405^2 + standard_dev_405^2)
% 
% % Go through each data point after index 35, going down and if it plus its
% % error bar is out of the threshold range, then that point is the knee











% 365 analysis

%%% FORWARD
% Calculate average of three runs
firstrun365 = forward_current_A365(1:301);
secondrun365 = forward_current_A365(302:602);
thirdrun365 = forward_current_A365(603:903);

averagecurrent365 = [];
current_standard_dev_365 = [];
for i = 1:301;
    averagecurrent365i = (firstrun365(i) + secondrun365(i) + thirdrun365(i))/3;
    current_standard_dev_365(i) = std([firstrun365(i); secondrun365(i); thirdrun365(i)]);
    averagecurrent365(i) =  averagecurrent365i;
end
voltagetouse365 = forward_voltage_V365(1:301);

% Uncertainty in Voltage
voltage_unc_365 = 0.001*ones(size(voltagetouse365));

% Uncertainty in Current
current_standard_error_365 = current_standard_dev_365/ sqrt(3);


% Find imax
figure
fitop = fitoptions('Method','NonlinearLeastSquares',...
               'StartPoint',[4e-8 4e-10 0.04]);
g = fittype('a-b*exp(-c*x)','options',fitop);
% g = fittype('a-b*exp(-c*x)');
f0 = fit(voltagetouse365,averagecurrent365',g);
f0_coeff = coeffvalues(f0);
plot(f0,voltagetouse365,averagecurrent365');






%%% REVERSE
% Calculate average of three runs 
firstrun_reverse_365 = reverse_current_A365(1:41);
secondrun_reverse_365 = reverse_current_A365(42:82);
thirdrun_reverse_365 = reverse_current_A365(83:123);

averagecurrent_reverse_365 =[];
current_reverse_standard_dev_365 = [];
for i = 1:41;
    averagecurrent_reverse_365i = (firstrun_reverse_365(i) + secondrun_reverse_365(i) + thirdrun_reverse_365(i)) / 3;
    current_reverse_standard_dev_365(i) = std([firstrun365(i); secondrun365(i); thirdrun365(i)]);
    averagecurrent_reverse_365(i) = averagecurrent_reverse_365i;
end
voltagetouse_reverse_365 = reverse_voltage_V365(1:41);
voltagetouse_reverse_365 = -voltagetouse_reverse_365;

% Uncertainty in Voltage
voltage_unc_reverse_365 = 0.001*ones(size(voltagetouse_reverse_365));

% Uncertainty in Current
current_reverse_standard_error_365 = current_reverse_standard_dev_365/ sqrt(3);




% Plot a figure with both data sets on one graph
figure
errorbar(voltagetouse365, averagecurrent365',current_standard_error_365,current_standard_error_365,voltage_unc_365,voltage_unc_365);
hold on;
errorbar(voltagetouse_reverse_365, averagecurrent_reverse_365,current_reverse_standard_error_365,current_reverse_standard_error_365,voltage_unc_reverse_365,voltage_unc_reverse_365)
plot(voltagetouse365, averagecurrent365');
plot(voltagetouse_reverse_365, averagecurrent_reverse_365);
title("Current vs. Voltage (365nm Filter)")
xlabel('Voltage (V)') 
ylabel('Current (A)')


% Plot just the forward bias
figure
errorbar(voltagetouse365, averagecurrent365',current_standard_error_365,current_standard_error_365,voltage_unc_365,voltage_unc_365);
title("Forward Bias Current vs. Voltage (365nm Filter)")
xlabel('Voltage (V)') 
ylabel('Current (A)')
hold on;


% Plot just the reverse bias 
figure
errorbar(voltagetouse_reverse_365, averagecurrent_reverse_365,current_reverse_standard_error_365,current_reverse_standard_error_365,voltage_unc_reverse_365,voltage_unc_reverse_365)
title("Reverse Bias Current vs. Voltage (365nm Filter)")
xlabel('Voltage (V)') 
ylabel('Current (A)')
hold on;


% Find knee using plateau uncertaities and deviation from that value

% Calculate aveage data point
running_total_365 = 0;
h365 = averagecurrent_reverse_365(35:40);
for i = 1:6;
    running_total_365 = h365(i) + running_total_365;
end
average_threshold_point_365 = running_total_365 / 6;

% Calculate average error bar size
running_total_365 = 0;
g365 = current_reverse_standard_error_365(35:40);
for i = 1:6;
    running_total_365 = g365(i)^2 + running_total_365;
end
uncertainty_mean_365 = sqrt(running_total_365)/6;

% Calculate standard deviation
standard_dev_365 = std(current_reverse_standard_error_365(35:40));

% Calculate threshold bar
threshold_bar_365 = sqrt(uncertainty_mean_365^2 + standard_dev_365^2);

% Calculate threshold
threshold = threshold_bar_365 + average_threshold_point_365;

% Plot the threshold line
figure
errorbar(voltagetouse_reverse_365, averagecurrent_reverse_365,current_reverse_standard_error_365,current_reverse_standard_error_365,voltage_unc_reverse_365,voltage_unc_reverse_365)
title("Reverse Bias Current vs. Voltage (365nm Filter)")
xlabel('Voltage (V)') 
ylabel('Current (A)')
hold on;
x = [-2 : 0.5 : 0];
Z = threshold * ones(1, length(x));
plot(x, Z, 'r', 'LineWidth', 2)

% Find knee using intersection of flat slopes

% Calculating slope of top line
k365 = averagecurrent_reverse_365(5:20);
j365 = voltagetouse_reverse_365(5:20);

% Plot intersection of lines
figure
errorbar(voltagetouse_reverse_365, averagecurrent_reverse_365,current_reverse_standard_error_365,current_reverse_standard_error_365,voltage_unc_reverse_365,voltage_unc_reverse_365)
title("Reverse Bias Current vs. Voltage (365nm Filter)")
xlabel('Voltage (V)') 
ylabel('Current (A)')
hold on;

% Plot top line
c = polyfit(j365,k365',1);
xFit = linspace(-2, 0, 100);
yFit = polyval(c, xFit);
hold on;
plot(xFit, yFit);

% Plot bottom line
x = [-2 : 0.5 : 0];
% Use already calculated value for average_threshold_point to get bottom
% line
Z = average_threshold_point_365 * ones(1, length(x));
plot(x, Z)
grid on;






% TODO
% Find i max
% Plot stopping voltage vs light frequency
% Chose colors of graphs