function [T] = stagwallthick(x,Ts,t,Ti,AblationRate)

k = 1.8; % thermal conductivity (W/(m*K)) -> 0.1441 W/(m*K) per Btu*in/(hr*ft^2*F)
rho = 265; % density (kg/m^3) -> 16 kg/m^3 per lbm/ft^3
cp = 1592;%200; % specific heat at constant pressure (J/(kg*K)) -> 4186.8 J/(kg*K) per Btu/(lbm*F)
a = k/(rho*cp); % thermal diffusivity
n = length(t);
T = zeros(n,1);
AbDist = zeros(n,1);
AbDist(1) = 0;
T(1) = Ti;

for ii = 2:n
    
    AbDist(ii) = AbDist(ii-1) + (AblationRate(ii) + AblationRate(ii-1))/2*(t(ii) - t(ii-1));
    ii = ii+1;
end

for i = 2:n 
    
    T(i) = erf((x-AbDist(i))/(2*sqrt(a*t(i))))*(Ti-Ts(i))+Ts(i);
    i = i+1;
end  
% plot(t,T,'b',t,Tbond*ones(length(t),1),'k--',t,SF*Tbond*ones(length(t),1),'r--')
% legend('Bond Temp','Bond Temp Limit','w/S.F.')



