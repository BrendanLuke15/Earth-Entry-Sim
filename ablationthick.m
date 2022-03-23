function [x_AB,AblationRate] = ablationthick(q_tot,Pt,t)

Data = [40	4.4	200	6.42;
        42	1.8	200	3.73;
        45	4.9	240	9.45;
        65	1.6	200	6.955;
        73	8.9	120	7.38;
        73	13.3 120 8.61;
        107	2.3	55 2.3;
        114	20.9 80 8.685;
        133	31.6 80 10.07;
        143	3.8	200	12.66;
        143	3.8	400	24.72;
        154	13.3 33	2.925;
        154	13.3 66	5.61;
        165	6.8	50 3.96;
        169	5 33 2.235;
        169	5 66 4.28;
        175	13.9 50 5.22;
        183	26.7 50 6.765];
% Heat Flux(X);  Stag. Pressure(Y);  Time;  Ablation
% W/cm^2;             kPa;            s;      mm

abRate = Data(:,4)./Data(:,3); % (mm/s)

options = fitoptions('poly11');
options.Lower = [0 0.0002 0.002];
options.Upper = [0 0.0004 0.004];
f = fit([Data(:,1), Data(:,2)], abRate,'poly11',options);

AblationRate = zeros(length(t),1);
for i = 1:length(t)
    AblationRate(i) = f(q_tot(i),Pt(i)/1000)/1000; % (m/s)
    i = i+1;
end

x_AB = trapz(t,AblationRate);
% (m)
end