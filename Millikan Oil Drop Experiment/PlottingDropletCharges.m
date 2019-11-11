drop_number = 1:13;
charge = [0.2427,
0.1271,
0.1338,
0.2682,
0.2693,
0.1761,
0.1104,
0.2241,
0.64252,
0.72425,
0.34792,
0.20199,
0.42052];
unc_charge = [0.0027548,
0.001422,
0.0015538,
0.004052,
0.004061,
0.005598,
0.0035085,
0.00713,
0.020467,
0.023062,
0.011073,
0.006429,
0.03389];



x1 = drop_number;
y1 = charge;
%dx1 = errangle;
dy1 = unc_charge;

%Plot data with x and y error bars
errorbar(x1,y1,dy1,'bo')

grid on

xlabel('oil drop number')
ylabel('Charge x 10^-18 (C)')

spacing = .16;
yticks(spacing*(1:100))

axis([1 14 0 1])

hold off