function [T0max] = stagwallthickroot(x,Ts,t,Ti,AblationRate)

k = 1.6; % thermal conductivity (W/(m*K))
rho = 265; % density (kg/m^3)
cp = 1592; % specific heat at constant pressure (J/(kg*K))
a = k/(rho*cp);
n = length(t);
T0 = zeros(n,1);
TBond = 250+273.15;


for ii = 2:n
    AbDist(1) = 0;
    AbDist(ii) = AbDist(ii-1) + (AblationRate(ii) + AblationRate(ii-1))/2*(t(ii) - t(ii-1));
    ii = ii+1;
end

for i = 2:n 
    T0(1) = Ti;
    T0(i) = erf((x-AbDist(i))/(2*sqrt(a*t(i))))*(Ti-Ts(i))+Ts(i)-TBond;
    i = i+1;
end
T0max = max(T0);




