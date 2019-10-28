wavelengths = [500,450,400,350,300,250]
stopping_potentials = [0.19,0.48,0.83,1.28,1.89,2.74]
frequencies = [];
i = 1;
for a = wavelengths
    disp(a)
    freq = (3*10^8)/(a * 10^-9);
    frequencies(i) = freq;
    i = i + 1;
end
hold on
plot(frequencies, stopping_potentials)
xlim([0 12*10^14])
ylim([-5 3])
disp(frequencies)
i = 1;
for b = frequencies -1
    run = frequencies (i+1) - frequencies(i)
    rise = stopping_potentials(i+1) - stopping_potentials(i)
end
% rise and run are always the same
x = -1:2*10^14:12*10^14
y = (rise/run)*x - 2.45
plot(x,y)