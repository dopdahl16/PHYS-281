TITLE= 'Select the file with the data you want to bring into MATLAB';
[filename,filepath] = uigetfile('*.*', TITLE); %Prompts the user to select a data file
full_filename = fullfile( filepath, filename );

%T = readtable(full_filename,'ReadVariableNames',true);

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





firstrun577 = forward_current_A577(1:301);
secondrun577 = forward_current_A577(302:602);
thirdrun577 = forward_current_A577(603:903);

averagecurrent577 = [];
for i = 1:301;
    averagecurrent577i = (firstrun577(i) + secondrun577(i) + thirdrun577(i))/3;
    averagecurrent577(i) =  averagecurrent577i;
end
voltagetouse577 = forward_voltage_V577(1:301);

figure
plot(voltagetouse577, averagecurrent577);
hold on;

p = polyfit(voltagetouse577, averagecurrent577', 2);
y1 = polyval(p,voltagetouse577);
plot(voltagetouse577, y1);
y2 = polyval(p,linspace(0,60,602));
plot(linspace(0,60,602),y2);
q = polyder(p);
r = roots(q);
hold off;



firstrun546 = forward_current_A546(1:301);
secondrun546 = forward_current_A546(302:602);
thirdrun546 = forward_current_A546(603:903);

averagecurrent546 = [];
for i = 1:301;
    averagecurrent546i = (firstrun546(i) + secondrun546(i) + thirdrun546(i))/3;
    averagecurrent546(i) =  averagecurrent546i;
end
voltagetouse546 = forward_voltage_V546(1:301);

figure
plot(voltagetouse546, averagecurrent546);
hold on;

p = polyfit(voltagetouse546, averagecurrent546', 2);
y1 = polyval(p,voltagetouse546);
plot(voltagetouse546, y1);
y2 = polyval(p,linspace(0,60,602));
plot(linspace(0,60,602),y2);
q = polyder(p);
r = roots(q);
hold off;



firstrun436 = forward_current_A436(1:301);
secondrun436 = forward_current_A436(302:602);
thirdrun436 = forward_current_A436(603:903);

averagecurrent436 = [];
for i = 1:301;
    averagecurrent436i = (firstrun436(i) + secondrun436(i) + thirdrun436(i))/3;
    averagecurrent436(i) =  averagecurrent436i;
end
voltagetouse436 = forward_voltage_V436(1:301);

figure
plot(voltagetouse436, averagecurrent436);
hold on;

p = polyfit(voltagetouse436, averagecurrent436', 2);
y1 = polyval(p,voltagetouse436);
plot(voltagetouse436, y1);
y2 = polyval(p,linspace(0,60,602));
plot(linspace(0,60,602),y2);
q = polyder(p);
r = roots(q);
hold off;



firstrun405 = forward_current_A405(1:301);
secondrun405 = forward_current_A405(302:602);
thirdrun405 = forward_current_A405(603:903);

averagecurrent405 = [];
for i = 1:301;
    averagecurrent405i = (firstrun405(i) + secondrun405(i) + thirdrun405(i))/3;
    averagecurrent405(i) =  averagecurrent405i;
end
voltagetouse405 = forward_voltage_V405(1:301);

figure
plot(voltagetouse405, averagecurrent405);
hold on;

p = polyfit(voltagetouse405, averagecurrent405', 2);
y1 = polyval(p,voltagetouse405);
plot(voltagetouse405, y1);
y2 = polyval(p,linspace(0,60,602));
plot(linspace(0,60,602),y2);
q = polyder(p);
r = roots(q);
hold off;



firstrun365 = forward_current_A365(1:301);
secondrun365 = forward_current_A365(302:602);
thirdrun365 = forward_current_A365(603:903);

averagecurrent365 = [];
for i = 1:301;
    averagecurrent365i = (firstrun365(i) + secondrun365(i) + thirdrun365(i))/3;
    averagecurrent365(i) =  averagecurrent365i;
end
voltagetouse365 = forward_voltage_V365(1:301);

figure
plot(voltagetouse365, averagecurrent365');
hold on;

p = polyfit(voltagetouse365, averagecurrent365', 2);
y1 = polyval(p,voltagetouse365);
plot(voltagetouse365, y1);
y2 = polyval(p,linspace(0,60,602));
plot(linspace(0,60,602),y2);
q = polyder(p);
r = roots(q);
hold off;





firstrun_reverse_577 = reverse_current_A577(1:41);
secondrun_reverse_577 = reverse_current_A577(42:82);
thirdrun_reverse_577 = reverse_current_A577(83:123);

averagecurrent_reverse_577 =[];
for i = 1:41;
    averagecurrent_reverse_577i = (firstrun_reverse_577(i) + secondrun_reverse_577(i) + thirdrun_reverse_577(i));
    averagecurrent_reverse_577(i) = averagecurrent_reverse_577i;
end
voltagetouse_reverse_577 = reverse_voltage_V577(1:41);
voltagetouse_reverse_577 = -voltagetouse_reverse_577;

figure
plot(voltagetouse_reverse_577, averagecurrent_reverse_577)
hold on;
p = polyfit(voltagetouse_reverse_577, averagecurrent_reverse_577', 25);
y1 = polyval(p,voltagetouse_reverse_577);
plot(voltagetouse_reverse_577, y1);



firstrun_reverse_546 = reverse_current_A546(1:41);
secondrun_reverse_546 = reverse_current_A546(42:82);
thirdrun_reverse_546 = reverse_current_A546(83:123);

averagecurrent_reverse_546 =[];
for i = 1:41;
    averagecurrent_reverse_546i = (firstrun_reverse_546(i) + secondrun_reverse_546(i) + thirdrun_reverse_546(i));
    averagecurrent_reverse_546(i) = averagecurrent_reverse_546i;
end
voltagetouse_reverse_546 = reverse_voltage_V546(1:41);
voltagetouse_reverse_546 = -voltagetouse_reverse_546;

figure
plot(voltagetouse_reverse_546, averagecurrent_reverse_546)
hold on;
p = polyfit(voltagetouse_reverse_546, averagecurrent_reverse_546', 25);
y1 = polyval(p,voltagetouse_reverse_546);
plot(voltagetouse_reverse_546, y1);



firstrun_reverse_436 = reverse_current_A436(1:41);
secondrun_reverse_436 = reverse_current_A436(42:82);
thirdrun_reverse_436 = reverse_current_A436(83:123);

averagecurrent_reverse_436 =[];
for i = 1:41;
    averagecurrent_reverse_436i = (firstrun_reverse_436(i) + secondrun_reverse_436(i) + thirdrun_reverse_436(i));
    averagecurrent_reverse_436(i) = averagecurrent_reverse_436i;
end
voltagetouse_reverse_436 = reverse_voltage_V436(1:41);
voltagetouse_reverse_436 = -voltagetouse_reverse_436;

figure
plot(voltagetouse_reverse_436, averagecurrent_reverse_436)
hold on;
p = polyfit(voltagetouse_reverse_436, averagecurrent_reverse_436', 25);
y1 = polyval(p,voltagetouse_reverse_436);
plot(voltagetouse_reverse_436, y1);



firstrun_reverse_405 = reverse_current_A405(1:41);
secondrun_reverse_405 = reverse_current_A405(42:82);
thirdrun_reverse_405 = reverse_current_A405(83:123);

averagecurrent_reverse_405 =[];
for i = 1:41;
    averagecurrent_reverse_405i = (firstrun_reverse_405(i) + secondrun_reverse_405(i) + thirdrun_reverse_405(i));
    averagecurrent_reverse_405(i) = averagecurrent_reverse_405i;
end
voltagetouse_reverse_405 = reverse_voltage_V405(1:41);
voltagetouse_reverse_405 = -voltagetouse_reverse_405;

figure
plot(voltagetouse_reverse_405, averagecurrent_reverse_405)
hold on;
p = polyfit(voltagetouse_reverse_405, averagecurrent_reverse_405', 25);
y1 = polyval(p,voltagetouse_reverse_405);
plot(voltagetouse_reverse_405, y1);



firstrun_reverse_365 = reverse_current_A365(1:41);
secondrun_reverse_365 = reverse_current_A365(42:82);
thirdrun_reverse_365 = reverse_current_A365(83:123);

averagecurrent_reverse_365 =[];
for i = 1:41;
    averagecurrent_reverse_365i = (firstrun_reverse_365(i) + secondrun_reverse_365(i) + thirdrun_reverse_365(i));
    averagecurrent_reverse_365(i) = averagecurrent_reverse_365i;
end
voltagetouse_reverse_365 = reverse_voltage_V365(1:41);
voltagetouse_reverse_365 = -voltagetouse_reverse_365;

figure
plot(voltagetouse_reverse_365, averagecurrent_reverse_365)
hold on;
p = polyfit(voltagetouse_reverse_365, averagecurrent_reverse_365', 25);
y1 = polyval(p,voltagetouse_reverse_365);
plot(voltagetouse_reverse_365, y1);