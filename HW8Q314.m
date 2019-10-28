%% HW8Q314.m
%
% Author: Daniel Opdahl
% Last modified: 10/28/19

x_axis_frequencies = -1:2*10^14:12*10^14;
wavelengths = [500,450,400,350,300,250];
stopping_potentials = [0.19,0.48,0.83,1.28,1.89,2.74];
frequencies = [];
i = 1;
for a = wavelengths
    freq = (3*10^8)/(a * 10^-9);
    frequencies(i) = freq;
    i = i + 1;
end

x1 = frequencies;
y1 = stopping_potentials;

%Curve fitting
myfittype1 = fittype('poly1');
myfit1 = fit(x1',y1',myfittype1);

coeffvals = coeffvalues(myfit1);
intercept = coeffvals(2);
slope = coeffvals(1);

work_function = -1*intercept % Units of eV
e = 1.60217662 * 10^-19;

i = 1;
tot = 0;
for voltage = stopping_potentials
    h = ((e* voltage) + work_function)/(frequencies(i));
    tot = tot + h;
    i = i + 1;
end

experimental_val_of_h = tot/6 % Units of eV*s