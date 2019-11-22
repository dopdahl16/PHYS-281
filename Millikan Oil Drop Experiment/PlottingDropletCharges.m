drop_number = 1:8;
charge = 10 * [
%     Insert data for charges here
];

unc_charge = 10^-1 * [
%     Insert data for charge uncertainty here
];

x1 = drop_number;
y1 = charge;
dy1 = unc_charge;

%Plot data with x and y error bars
errorbar(x1,y1,dy1,'bo')

grid on

xlabel('Oil Drops')
ylabel('Charge x 10^-19 (C)')

spacing = .15;
yticks(spacing*(1:100))

axis([0 9 0 8])

hold off