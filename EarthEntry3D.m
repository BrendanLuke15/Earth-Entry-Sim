% Brendan Luke
% November 2, 2021
clear, clc
format compact
close all
tic
%% Input Parameters
% input values: CONSTANT
R0 = 6378137; % radius of Earth (m)
G = 6.6743*10^-11; % gravitational constant (N*m^2/kg^2)
M = 5.97237*10^24; % mass of Earth (kg)
wEarth = 72.92115*10^-6; % angular velocity of Earth rotation (rad/s)
f = 1/298.257223563; % WGS-84 flattenig ellipsoid
EarthImg = imread('eo_base_2020_clean_3600x1800.png');
% EarthImg = imread('eo_base_2020_clean_720x360.jpg');
EarthImg = flip(EarthImg,1); % rectify image

% input values: DESIGN
m = 111; % mass of entry vehicle (kg)
Cd = 2.1; % drag coefficient
Cl = 0; % lift coefficient
d = 0.625; % entry vehicle diameter (m)
S = (d/2)^2*pi(); % ref. area (m^2)
BC = m/(S*Cd); % ballistic coefficient (kg/m^2)
Rn = 0.22; % nose radius (m)
emis = 0.9; % thermal IR emmisivity (0.9 for PICA)
Conditions = [22,418]; % periapsis and apoapsis altitudes (km)

% orbit fit data
%{
dt	h	f(t)	Abs error %		Amplitude	Frequency	Phase
8.55	331.454	331.6413089	0.1%		198.3483836	0.070678109	0.369382132
12.78333333	279.966	278.3381586	0.6%				
17.55	212.388	212.3869411	0.0%		Error Sum:		
20.3	173.772	174.2519538	0.3%		0.9%		
23.51666667	131.938	131.937254	0.0%				

%}

% Initial Orbit Angles
inc = 51.6; % inclination (°)
RAAN = 251; % right ascension of ascending node (°)
w = 64; % argument of periapsis (°)

fprintf('The entry vehicle has a ballistic coefficient of %.1f kg/m^2.\n',BC)

%% Atmosphere Model
h = 0:1:120000;
for i = 1:length(h)
    A = EarthAtmos(h(i),R0);
    T(i) = A(2);
    P(i) = A(1);
    rho(i) = A(3);
end
figure('name','Atmosphere Model');
set(gcf,'WindowState','maximized');
yyaxis right
plot(h/1000,T,'LineWidth',2);
ylabel('Temperature');
yyaxis left
semilogy(h/1000,rho,'LineWidth',2);
hold on
semilogy(h/1000,P/101325,'-','Color',[0.9290, 0.6940, 0.1250],'LineWidth',2);
ylabel('Density and Pressure');
xlabel('Altitude (km)');
legend('Density (kg/m^3)','Pressure (atm)','Temperature (K)');
grid on

%% Orbit Entry Paramters
a = 1000*mean(Conditions)+R0; % orbit semi major axis (m)
V_Entry = sqrt(G*M*(2/(120000+R0)-1/a)); % entry velocity (m/s)
angMomentum = (Conditions(1)*1000+R0)*sqrt(G*M*(2/(Conditions(1)*1000+R0)-1/a)); % angular momentum m^2/s
p = angMomentum^2/(G*M); % semi-latus rectum (m)
e = (1000*diff(Conditions))/(2*R0+1000*sum(Conditions)); % orbital eccentricity
EFPA = -acosd(angMomentum/V_Entry/(120000+R0)); % entry flight path angle (°)
fprintf('The entry velocity is %.0f m/s and the EFPA is %.2f° (Inertial).\n',V_Entry,EFPA)

% Rotation Matrix
R11 = cosd(RAAN)*cosd(w) - sind(RAAN)*sind(w)*cosd(inc); 
R12 = -cosd(RAAN)*sind(w) - sind(RAAN)*cosd(w)*cosd(inc);
R13 = sind(RAAN)*sind(i);
R21 = sind(RAAN)*cosd(w) + cosd(RAAN)*sind(w)*cosd(inc);
R22 = -sind(RAAN)*sind(w) + cosd(RAAN)*cosd(w)*cosd(inc);
R23 = -cosd(RAAN)*sind(inc);
R31 = sind(w)*sind(inc);
R32 = cosd(w)*sind(inc);
R33 = cosd(inc);

R = [R11 R12 R13;
     R21 R22 R23;
     R31 R32 R33];

% create data for plot
for i = 1:3600
    r_plot(i) = p/(1+e*cosd(i/10));
    x_plot = r_plot(i)*cosd(i/10);
    y_plot = r_plot(i)*sind(i/10);
    z_plot = 0;
    
    temp = (R*[x_plot; y_plot; z_plot])';
    xSetup(i) = temp(1);
    ySetup(i) = temp(2);
    zSetup(i) = temp(3);    
end

% Get Orbit details
[~,idx_max] = max(r_plot); % apoapsis
[~,idx_min] = min(r_plot); % periapsis
[~,idx_Entry] = min(abs(r_plot(1800:end)-(R0+120000)));
idx_Entry = idx_Entry+1800; % descending half of orbit
EccenAnom = asind((sqrt(1-e^2)*sind(idx_Entry/10))/(1+e*cosd(idx_Entry/10)));
V_Entry_vect = R*[-V_Entry*sind(EccenAnom); V_Entry*(sqrt(1-e^2)*cosd(EccenAnom)); 0]; % vector entry velocity (inertial)

% make plot
figure('name','Entry Geometry');
set(gcf,'WindowState','maximized');
[Xs,Ys,Zs] = ellipsoid(0,0,0,R0,R0,R0*(1-f),360);
surf(Xs,Ys,Zs);
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',EarthImg,'EdgeColor',[0 0 0],'EdgeAlpha',0);%.3);
hold on
plot3(xSetup,ySetup,zSetup,'LineWidth',2);
plot3(xSetup(idx_max),ySetup(idx_max),zSetup(idx_max),'k.','MarkerSize',20); % apoapsis point
plot3(xSetup(idx_min),ySetup(idx_min),zSetup(idx_min),'k.','MarkerSize',20); % periapsis point
text(1.08*xSetup(idx_max),1.08*ySetup(idx_max),1.08*zSetup(idx_max),['Apoapsis ',num2str(Conditions(2)),' km'],'HorizontalAlignment','Right');
text(1.08*xSetup(idx_min),1.08*ySetup(idx_min),1.08*zSetup(idx_min),['Periapsis ',num2str(Conditions(1)),' km']);
plot3(xSetup(idx_Entry),ySetup(idx_Entry),zSetup(idx_Entry),'g.','MarkerSize',30); % entry interface point
text(1.08*xSetup(idx_Entry),1.08*ySetup(idx_Entry),1.08*zSetup(idx_Entry),'Entry Interface (120 km)');
plot3([xSetup(idx_Entry) xSetup(idx_Entry)+100*V_Entry_vect(1)],[ySetup(idx_Entry) ySetup(idx_Entry)+100*V_Entry_vect(2)],[zSetup(idx_Entry) zSetup(idx_Entry)+100*V_Entry_vect(3)],'r','LineWidth',2); % entry velocity vector
grid on
axis equal
axis off
set(gca,'clipping','off');

%% Trajectory Simulation
z0 = [xSetup(idx_Entry); ySetup(idx_Entry); zSetup(idx_Entry); V_Entry_vect(1); V_Entry_vect(2); V_Entry_vect(3)];

ReEntryEventsFcn2 = @(t,z) ReEntryEventsFcn(t,z,R0,f);
options = odeset('MaxStep',1,'Stats','off','RelTol',1e-9,'Events',ReEntryEventsFcn2);
[t,z] = ode45(@(t,z) dz(t,z,G,M,BC,R0,Cl,Cd,wEarth,f),[0 1800],z0,options);

%% Process Data
for i = 1:length(t)
    temp2 = dz(t(i),z(i,:),G,M,BC,R0,Cl,Cd,wEarth,f); % slopes
    [longitude(i),latitude_geo(i),radius(i)] = cart2sph(z(i,1), z(i,2), z(i,3));
    longitude_Plot(i) = longitude(i) - t(i)*wEarth; % correct Longitude for rotating Earth
    velocity_Inertial(i) = norm([z(i,4), z(i,5), z(i,6)]);
    altitude(i) = radius(i) - R0^3*(1-f)/sqrt(R0^4*(1-f)^2*(cos(longitude(i))^2+sin(longitude(i))^2)*cos(latitude_geo(i))^2+R0^4*sin(latitude_geo(i))^2);
    velocity_Air(i) = norm([z(i,4), z(i,5), z(i,6)] - wEarth*radius(i)*[-sin(longitude(i))*cos(latitude_geo(i)), cos(longitude(i))*cos(latitude_geo(i)), 0]); % air speed
    temp = EarthAtmos(altitude(i),R0); % atmospherics
    aG(i) = G*M/(radius(i)^2);
    aD(i) = EarthAtmosRho(altitude(i),R0)/2/BC*velocity_Air(i)^2;
    aL(i) = EarthAtmosRho(altitude(i),R0)/2/BC/Cd*Cl*velocity_Air(i)^2;
    Q(i) = 1/2*temp(3)*velocity_Air(i)^2;
    Pt(i) = Q(i) + temp(1);
    Mach(i) = velocity_Air(i)/temp(5);
    qConv(i) = 1.7415*10^-4*sqrt(temp(3)/Rn)*velocity_Air(i)^3/100/100;
    qRad(i) = 8*2.787*10^-67*Rn^0.2*temp(3)^1.05*velocity_Air(i)^19/100/100;
    qTot(i) = qConv(i) + qRad(i);
    %
    angMo = norm(cross([z(i,1), z(i,2), z(i,3)],[z(i,4), z(i,5), z(i,6)] - wEarth*radius(i)*[-sin(longitude(i))*cos(latitude_geo(i)), cos(longitude(i))*cos(latitude_geo(i)), 0])); % angular momentum m^2/s
    FPA(i) = -acosd(angMo/velocity_Inertial(i)/radius(i)); % flight path angle, not inertial
    vertAccel(i) = aD(i)*sind(-FPA(i)) + aL(i)*cosd(FPA(i));%dot([temp2(4), temp2(5), temp2(6)],[z(i,1), z(i,2), z(i,3)])/radius(i);
    vertVelo_Inertial(i) = dot([z(i,4), z(i,5), z(i,6)],[z(i,1), z(i,2), z(i,3)])/radius(i);
    %}
    GForce(i) = sqrt(aD(i)^2+aL(i)^2)/9.81;
    if i == 1
        HeatLoad(i) = 0;
        angle(i) = 0;
        xPos(i) = 0;
        yPos(i) = z(i,2);
        Ti = temp(2);
    else
        HeatLoad(i) = HeatLoad(i-1) + (qTot(i) + qTot(i-1))/2*(t(i) - t(i-1));
        angle(i) = angle(i-1) + 2*(z(i,1)-z(i-1,1))/(z(i,2)+z(i-1,2));
        xPos(i) = z(i,2)*sin(angle(i));
        yPos(i) = z(i,2)*cos(angle(i));
    end
    Tw(i) = ((100^2*qTot(i))/(5.67*10^-8*emis))^0.25 + temp(2);
%     Tw(i) = ((100^2*qTot(i))/(5.67*10^-8*emis) + temp(2)^4)^0.25;
    Tsur(i) = temp(2);  
end

% ablation: recession distance
[x_AB,AblationRate] = ablationthick(qTot,Q,t); % (m)

% Insulation
Tbond = (250+273.15);%644; % (K) from holy grail slideshow
x_tot = round(x_AB,3); % (m)
Ti = ((0.85*230)/((5.67*10^-8*emis)))^0.25; % (k) inital temp of material, calculated from re-radiation temp. from Earth IR with emissivity (0.9) and reflectivity (0.85)

while 1
    T = stagwallthick(x_tot,Tw,t,Ti,AblationRate);
    if max(T) > Tbond
        x_tot = x_tot+0.001; % increment millimeter at a time
    else
        break
    end
end

% x_tot = x_tot-round(x_AB,3);
T_bond = stagwallthick(x_tot,Tw,t,Ti,AblationRate);

figure('name','Bondline Temperature');
set(gcf,'WindowState','maximized');
plot(t,T_bond,'b','LineWidth',1.5)
hold on
plot(t,Tbond*ones(length(t),1),'r--','LineWidth',1.5)
ylabel('Temperature (K)');
xlabel('Time (s)');
title(['Nominal Thickness: ',num2str(100*x_tot),' cm']);
legend('Bond Temp.','Bond Temp. Limit');%,'Location','East');
xlim([0 t(end)]);
grid on

PICA_Envelope = [0 1500;      % known envelope
                 1 1500;
                 1 0;
                 0 1800;      % newly shown envelope
                 1.283 1800;
                 1.283 0]; 
             
fprintf('The maximum inertial loading is %.1f gs.\n',max(GForce))
fprintf('The peak heating is %.1f W/cm^2.\n',max(qTot))
fprintf('The total heat load is %.0f J/cm^2.\n',max(HeatLoad))
fprintf('Zero-margin ablated thickness: %.1f cm.\n',100*x_AB)
fprintf('Zero-margin total thickness: %.1f cm.\n',100*x_tot)

landingLat_geo = 180/pi*latitude_geo(end);
landingLong = 180/pi*(longitude(end)-wEarth*t(end));

fprintf('Landing Site: %0.1f° N, %0.1f° E.\n',landingLat_geo,landingLong)

%% Plot Data
figure('name','TPS Performance Map');
set(gcf,'WindowState','maximized');
plot(Q/101325,qTot,'LineWidth',2);
hold on
plot(PICA_Envelope(1:3,1),PICA_Envelope(1:3,2),'k--','HandleVisibility','off','LineWidth',1.5);
plot(PICA_Envelope(4:6,1),PICA_Envelope(4:6,2),'k--','HandleVisibility','off','LineWidth',1.5);
text(0.1,1400,'PICA Envelope');
text(0.3,1650,'PICA Envelope Expansion');
axis([0 1.4 0 2000]);
xlabel('Stag. Pressure (atm)');
ylabel('Stag. Heat Flux (W/cm^{2})');
grid on

figure('name','G-Force');
set(gcf,'WindowState','maximized');
plot(t,GForce,'LineWidth',2);
grid on
xlabel('Time (s)');
ylabel("G-Froce (g's)");

%{
figure('name','Vertical Profile');
set(gcf,'WindowState','maximized');
yyaxis left
plot(z(:,1)/1000,z(:,2)/1000-R0/1000,'LineWidth',2);
ylabel('Altitude (km)');
yyaxis right
plot(z(:,1)/1000,z(:,3).*sin(z(:,4)),'LineWidth',2);
ylabel('Vertical Velocity (m/s)');
xlabel('Downrange Distance (km)');
grid on
%}

%{
figure('name','Vertical Acceleration');
set(gcf,'WindowState','maximized');
plot(t,vertAccel,'LineWidth',2);
xlabel('Time (s)');
ylabel('Vertical Acceleration (m/s^2)');
grid on
%}

figure('name','Loading Profile');
set(gcf,'WindowState','maximized');
plot(GForce,altitude/1000,'LineWidth',2);
xlabel('G-Force (g)');
ylabel('Altitude (km)');
ylim([0 120]);
grid on

figure('name','Heating');
set(gcf,'WindowState','maximized');
yyaxis right
plot(t,qConv,'g-','LineWidth',1.5)
hold on
plot(t,qRad,'r-','LineWidth',1.5)
plot(t,qTot,'b-','LineWidth',1.5)
ylabel('Heat Rate (W/cm^2)');
hold on
yyaxis left
plot(t,Tw,'k','LineWidth',1.5)
hold on
plot(t,Tsur,'--k','LineWidth',1.5)
% yline(Tmax,'-r','LineWidth',2.5);
% yline(1810,'--r','LineWidth',2.5);
ylabel('Temperature (K)');
xlabel('Time (s)');
title('Stagnation Point Thermal Environment');
legend('Radiative Equilibrium','Surrounding','Convective','Radiative','Total');
grid on

figure('name','Heat Load');
set(gcf,'WindowState','maximized');
yyaxis right
plot(t,HeatLoad,'k-','LineWidth',2)
ylabel('Heat Load (J/cm^2)');
hold on
yyaxis left
plot(t,qConv,'g-',t,qRad,'r-',t,qTot,'b-','LineWidth',2)
ylabel('Heat Rate (W/cm^2)');
xlabel('Time (s)');
legend('Convective','Radiative','Total','Heat Load','Location','East');
grid on

figure('name','Altitude VS Air Speed');
set(gcf,'WindowState','maximized');
plot(velocity_Air,altitude/1000,'LineWidth',2);
xlabel('Air Speed (m/s)');
ylabel('Altitude (km)');
ylim([0 120]);
grid on

figure('name','Simulated Trajectory');
set(gcf,'WindowState','maximized');
[Xs,Ys,Zs] = ellipsoid(0,0,0,R0,R0,R0*(1-f),360);
surf(Xs,Ys,Zs);
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',EarthImg,'EdgeColor',[0 0 0],'EdgeAlpha',0);%.3);
hold on
% plot3(z(:,1),z(:,2),z(:,3),'LineWidth',2); % trajectory
for i = 1:length(t)
    [x_traj(i),y_traj(i),z_traj(i)] = sph2cart(longitude_Plot(i),latitude_geo(i),radius(i));
end
plot3(x_traj,y_traj,z_traj,'LineWidth',2); % trajectory
plot3(xSetup(idx_Entry),ySetup(idx_Entry),zSetup(idx_Entry),'g.','MarkerSize',30); % entry interface point
text(1.08*xSetup(idx_Entry),1.08*ySetup(idx_Entry),1.08*zSetup(idx_Entry),'Entry Interface (120 km)');
[xL,yL,zL] = sph2cart(pi/180*landingLong,pi/180*landingLat_geo,radius(end));
plot3(xL,yL,zL,'r.','MarkerSize',20); % landing site point
textstr = sprintf('Landing Site: %0.1f° N, %0.1f° E.\n',landingLat_geo,landingLong);%fprintf('Landing Site: '+num2str(landingLat_geo,'%0.1f')+'° N, '+num2str(landingLong,'%0.1f')+'° E');
text(1.08*xL,1.08*yL,1.08*zL,textstr);
grid on
axis equal
axis off
set(gca,'clipping','off');

%% Fancy SE Plot
% find values
[~,maxMach] = max(Mach);
[~,maxG] = max(GForce);
[~,maxHeat] = max(qTot);
[~,Mach1] = min(abs(Mach-1));
[~,parachute] = min(abs(altitude-5500));
[~,vertVelo] = min(vertVelo_Inertial);
[~,maxVertAccel] = max(vertAccel);

% plot
figure('name','SE Plot','windowstate','maximized');
xlabel('Time, E+ (mins)');
grid on
xlim([0 14]);
% altitude side
yyaxis left
plot(t(1:parachute)/60,altitude(1:parachute)/1000,'LineWidth',2);
hold on
plot(t(maxG)/60,altitude(maxG)/1000,'k.','MarkerSIze',25);
text(t(maxG)/60-0.3,altitude(maxG)/1000-2,['E+'+string(minutes(t(maxG)/60),"mm:ss")+newline+'Max g-force: '+num2str(max(GForce),"%.1f")],...
    'horizontalalignment','center','verticalalignment','top');
plot(t(parachute)/60,altitude(parachute)/1000,'k.','MarkerSIze',25);
text(t(parachute)/60+0.8,altitude(parachute)/1000+2,['E+'+string(minutes(t(parachute)/60),"mm:ss")+newline+'Parachute Deploy (~5500 m)'],...
    'horizontalalignment','center','verticalalignment','bottom');
plot(t(vertVelo)/60,altitude(vertVelo)/1000,'k.','MarkerSIze',25);
text(t(vertVelo)/60+0.8,altitude(vertVelo)/1000+2,['E+'+string(minutes(t(vertVelo)/60),"mm:ss")+newline+'Max vertical speed: '+num2str(abs(min(vertVelo_Inertial)),"%.0f")+' m/s'],...
    'horizontalalignment','center','verticalalignment','bottom');
ylabel('Altitude (km)');

% velocity side
yyaxis right
plot(t(1:parachute)/60,velocity_Air(1:parachute),'LineWidth',2);
hold on
plot(t(maxMach)/60,velocity_Air(maxMach),'k.','MarkerSIze',25);
text(t(maxMach)/60,velocity_Air(maxMach)-50,['E+'+string(minutes(t(maxMach)/60),"mm:ss")+newline+'Max Mach: M='+num2str(max(Mach),"%.0f")],...
    'horizontalalignment','center','verticalalignment','top');
plot(t(Mach1)/60,velocity_Air(Mach1),'k.','MarkerSIze',25);
text(t(Mach1)/60+0.05,velocity_Air(Mach1)+50,['E+'+string(minutes(t(Mach1)/60),"mm:ss")+newline+'Mach 1'],...
    'horizontalalignment','center','verticalalignment','bottom');
plot(t(maxHeat)/60,velocity_Air(maxHeat),'k.','MarkerSIze',25);
text(t(maxHeat)/60-0.4,velocity_Air(maxHeat)-50,['E+'+string(minutes(t(maxHeat)/60),"mm:ss")+newline+'Max Heating: '+num2str(max(qTot),"%.0f")+' W/cm^2'],...
    'horizontalalignment','center','verticalalignment','top');
ylabel('Air Speed (m/s)');

%% Stop Clock
fprintf('\n');
toc

function dzdt = dz(t,z,G,M,BC,R0,Cl,Cd,wEarth,f)
% z1 = x, z2 = y, z3 = z,
% z4 = vx, z5 = vy, z6 = vz
[long,lat,r] = cart2sph(z(1), z(2), z(3));
v_Air = [z(4), z(5), z(6)] - wEarth*norm(r)*[-sin(long)*cos(lat), cos(long)*cos(lat), 0]; % air speed
altitude = r - R0^3*(1-f)/sqrt(R0^4*(1-f)^2*(cos(long)^2+sin(long)^2)*cos(lat)^2+R0^4*sin(lat)^2);

%v_Air = [z(4), z(5), z(6)];
%norm(v_Air)

aG = G*M/(norm(r)^2);
aD = EarthAtmosRho(altitude,R0)/2/BC*norm(v_Air)^2;
aL = aD/Cd*Cl; % lift

h = cross([z(1), z(2), z(3)],v_Air); % angular momentum vector
liftVect = cross(v_Air,h); % lift vector

dzdt(1) = z(4);
dzdt(2) = z(5);
dzdt(3) = z(6);
dzdt(4) = -aG*z(1)/norm(r) - aD*v_Air(1)/norm(v_Air) + aL*liftVect(1)/norm(liftVect); % gravity (r), drag (v), lift (norm to v)
dzdt(5) = -aG*z(2)/norm(r) - aD*v_Air(2)/norm(v_Air) + aL*liftVect(2)/norm(liftVect);
dzdt(6) = -aG*z(3)/norm(r) - aD*v_Air(3)/norm(v_Air) + aL*liftVect(3)/norm(liftVect);
dzdt = dzdt';
end

% Stop at Zero Altitude
function [position,isterminal,direction] = ReEntryEventsFcn(t,z,R0,f)
[long,lat,r] = cart2sph(z(1), z(2), z(3));
altitude = r - R0^3*(1-f)/sqrt(R0^4*(1-f)^2*(cos(long)^2+sin(long)^2)*cos(lat)^2+R0^4*sin(lat)^2);

position = [altitude-120000 altitude]; % The value that we want to be zero, altitude or skip out case
isterminal = [1 1];  % Halt integration 
direction = [1 0];   % The zero can be approached from either direction for altitude, increasing for skipout
end