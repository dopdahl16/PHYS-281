%This script automates the data collection capabilities for the
%PASCO Photoelectric effect apparatus.

%Both the Keithley Picoammeter and the AMREL powersupply should be
%connected in series to the computer using the Agilent GPIB-to-USB cable, and the
%GPIB-to-GPIB cable.

%% Code exported from tmtool for Keithley Ammeter
% Find a VISA-GPIB object.
obj1 = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::14::INSTR', 'Tag', '');

% Create the VISA-GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('AGILENT', 'GPIB0::14::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1);
end

% Connect to instrument object, obj1.
fopen(obj1);

%% Code exported from tmtool for AMREL powersupply
% Find a VISA-GPIB object.
obj2 = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB0::12::0::INSTR', 'Tag', '');

% Create the VISA-GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(obj2)
    obj2 = visa('AGILENT', 'GPIB0::12::0::INSTR');
else
    fclose(obj2);
    obj2 = obj2(1);
end

% Connect to instrument object, obj2.
fopen(obj2);

%% Automating Photoelectric effect data acquistion
%Prompts for the user
vstart = input('Input your desired initial value for the powersupply voltage (in V): ');
vend = input('Input your desired final value for the powersupply voltage(in V): ');
vstep = input('Input the step size for the powersupply voltage (in V): ');

n = (vend - vstart)/vstep + 1;
voltage = (vstart:vstep:vend)';

%initialize output data vector
current = zeros(n,1);

for i = 1:n
    % Communicating with instrument object, obj2: Setting voltage
    fprintf(obj2, horzcat('vset ', num2str(voltage(i))));
    pause(1.0); %pause for 1 sec for current to change
    % Communicating with instrument object, obj1: Reading current
    tempcurrent = fscanf(obj1);
    cutcurrent = tempcurrent(5:end);
    current(i) = str2double(cutcurrent);

end

data = [voltage, current];

fprintf(obj2, 'vset 0');


volt_list = data(:,1);
current_list = data(:,2);

grid on

hold off

plot(volt_list, current_list, 'o');