% Brendan Luke
% February 16, 2021
% 3D Earth Entry Locating
clear, clc, close all
format compact
tic
%% Input Parameters
% input values: CONSTANT
R0 = 6378137; % radius of Earth (m)
G = 6.6743*10^-11; % gravitational constant (N*m^2/kg^2)
M = 5.97237*10^24; % mass of Earth (kg)
wEarth = 72.92115*10^-6; % angular velocity of Earth rotation (rad/s)
f = 1/298.257223563; % WGS-84 flattening ellipsoid
EarthImg = imread('eo_base_2020_clean_3600x1800.png');
% EarthImg = imread('eo_base_2020_clean_720x360.jpg');
EarthImg = flip(EarthImg,1); % rectify image

% input values: DESIGN
m = 21200/2.205; % mass of entry vehicle (kg)
Cd = 1.23; % drag coefficient
Cl = 0.13*Cd;%0.48*0.27*Cd; % lift coefficient
d = 4;%13/3.28; % entry vehicle diameter (m)
S = (d/2)^2*pi(); % ref. area (m^2)
BC = m/(S*Cd); % ballistic coefficient (kg/m^2)
Rn = 2; % nose radius (m)
emis = 0.9; % thermal IR emmisivity (0.9 for PICA)
Conditions = [22,418]; % periapsis and apoapsis altitudes (km)

% Initial Orbit Angles
inc = 51.6; % inclination (°)
% South Track
% RAAN = 253; % right ascension of ascending node (°)
% w = 65; % argument of periapsis (°)
% North Track
RAAN = 127; % right ascension of ascending node (°)
w = 165; % argument of periapsis (°)
% West Coast N
% RAAN = 94; % right ascension of ascending node (°)
% w = 156; % argument of periapsis (°)
% West Coast S
% RAAN = 210; % right ascension of ascending node (°)
% w = 70; % argument of periapsis (°)

fprintf('The entry vehicle has a ballistic coefficient of %.1f kg/m^2.\n',BC)

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
[Xs,Ys,Zs] = ellipsoid(0,0,0,R0,R0,R0*(1-f),36);
surf(Xs,Ys,Zs);
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',EarthImg,'EdgeColor',[0 0 0],'EdgeAlpha',0.3);
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
options = odeset('MaxStep',10,'Stats','off','RelTol',1e-9,'Events',ReEntryEventsFcn2);
[t,z] = ode45(@(t,z) dz(t,z,G,M,BC,R0,Cl,Cd,wEarth,f),[0 1800],z0,options);

%% Process Data
options = optimset('Display','off');
for i = 1:length(t)
    temp2 = dz(t(i),z(i,:),G,M,BC,R0,Cl,Cd,wEarth,f); % slopes
    [longitude(i),latitude_geo(i),radius(i)] = cart2sph(z(i,1), z(i,2), z(i,3));
    velocity_Inertial(i) = norm([z(i,4), z(i,5), z(i,6)]);
    altitude(i) = radius(i) - R0^3*(1-f)/sqrt(R0^4*(1-f)^2*(cos(longitude(i))^2+sin(longitude(i))^2)*cos(latitude_geo(i))^2+R0^4*sin(latitude_geo(i))^2);
    velocity_Air(i) = norm([z(i,4), z(i,5), z(i,6)] - wEarth*radius(i)*[-sin(longitude(i))*cos(latitude_geo(i)), cos(longitude(i))*cos(latitude_geo(i)), 0]); % air speed
    temp = EarthAtmos(altitude(i),R0); % atmospherics
    aG(i) = G*M/(radius(i)^2);
    aD(i) = EarthAtmosRho(altitude(i),R0)/2/BC*velocity_Air(i)^2;
    aL(i) = aD(i)/Cd*Cl;
    Q(i) = 1/2*temp(3)*velocity_Air(i)^2;
    Pt(i) = Q(i) + temp(1);
    Mach(i) = velocity_Air(i)/temp(5);
    qConv(i) = 1.7415*10^-4*sqrt(temp(3)/Rn)*velocity_Air(i)^3/100/100;
    qRad(i) = 8*2.787*10^-67*Rn^0.2*temp(3)^1.05*velocity_Air(i)^19/100/100;
    qTot(i) = qConv(i) + qRad(i);
    %{
    vertAccel(i) = temp2(3)*sin(z(i,4)) + z(i,3)*temp2(4)*cos(z(i,4));
    horiAccel(i) = temp2(3)*cos(z(i,4)) - z(i,3)*temp2(4)*sin(z(i,4));
    %}
    GForce(i) = sqrt(aD(i)^2+aL(i)^2)/9.81;
    longitude_Plot(i) = longitude(i)-wEarth*t(i);
    [x,y,zp] = sph2cart(longitude_Plot(i),latitude_geo(i),radius(i));
    x_plot(i) = x;
    y_plot(i) = y;
    z_plot(i) = zp;
    viewDist(i) = real(fsolve(@(x) horizontalViewDist(x,R0,altitude(i)), R0,options));
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
    %Tw(i) = ((100^2*qTot(i))/(5.67*10^-8*emis) + temp(2)^4)^0.25;
    Tw(i) = ((100^2*qTot(i))/(5.67*10^-8*emis))^0.25 + temp(2);
    Tsur(i) = temp(2);  
end

PICA_Envelope = [0 1500;      % known envelope
                 1 1500;
                 1 0;
                 0 1800;      % newly shown envelope
                 1.283 1800;
                 1.283 0]; 
             
fprintf('The maximum inertial loading is %.1f gs.\n',max(GForce))
fprintf('The peak heating is %.1f W/cm^2.\n',max(qTot))
fprintf('The total heat load is %.0f J/cm^2.\n',max(HeatLoad))

landingLat_geo = 180/pi*latitude_geo(end);
landingLong = 180/pi*(longitude(end)-wEarth*t(end));

fprintf('Landing Site: %0.1f° N, %0.1f° E.\n',landingLat_geo,landingLong)
fprintf('Max Wall Temp: %0.0f K.\n',max(Tw))

%% Plot Data

%{
figure('name','Heating');
% set(gcf,'WindowState','maximized');
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
%}

figure('name','Simulated Trajectory');
set(gcf,'WindowState','maximized');
[Xs,Ys,Zs] = ellipsoid(0,0,0,R0,R0,R0*(1-f),36);
surf(Xs,Ys,Zs);
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',EarthImg,'EdgeColor',[0 0 0],'EdgeAlpha',0.3);
hold on
plot3(x_plot,y_plot,z_plot,'LineWidth',2); % trajectory
plot3(xSetup(idx_Entry),ySetup(idx_Entry),zSetup(idx_Entry),'g.','MarkerSize',30); % entry interface point
text(1.08*xSetup(idx_Entry),1.08*ySetup(idx_Entry),1.08*zSetup(idx_Entry),'Entry Interface (120 km)');
[xL,yL,zL] = sph2cart(pi/180*landingLong,pi/180*landingLat_geo,radius(end));
plot3(xL,yL,zL,'r.','MarkerSize',20); % landing site point
text(1.08*xL,1.08*yL,1.08*zL,'Landing Site');
grid on
axis equal
axis off
set(gca,'clipping','off');

%
figure('name','Hot Zones');
set(gcf,'WindowState','maximized');
[Xs,Ys,Zs] = ellipsoid(0,0,0,R0,R0,R0*(1-f),36);
surf(Xs,Ys,Zs);
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',EarthImg,'EdgeColor',[0 0 0],'EdgeAlpha',0.3);
hold on
plot3(xSetup(idx_Entry),ySetup(idx_Entry),zSetup(idx_Entry),'g.','MarkerSize',30); % entry interface point
text(1.08*xSetup(idx_Entry),1.08*ySetup(idx_Entry),1.08*zSetup(idx_Entry),'Entry Interface (120 km)');
[xL,yL,zL] = sph2cart(pi/180*landingLong,pi/180*landingLat_geo,radius(end));
plot3(xL,yL,zL,'r.','MarkerSize',20); % landing site point
text(1.08*xL,1.08*yL,1.08*zL,'Landing Site');
colours = hot(100);
for i = 1:length(t)
    plot3(x_plot(i),y_plot(i),z_plot(i),'.','Color',colours(round(99*qTot(i)/max(qTot)+1),:),'MarkerSize',10); % trajectory points colour mapped to heat rate
end
grid on
axis equal
axis off
set(gca,'clipping','off');
%}

% 2D Ground Track Plot
figure('name','Ground Track');
set(gcf,'WindowState','maximized');
image([-180 180],[-90 90],flip(EarthImg,1));
hold on
for i = 1:length(t)
    plot(longitude_Plot(i)*180/pi,-latitude_geo(i)*180/pi,'.','Color',colours(round(99*qTot(i)/max(qTot)+1),:),'MarkerSize',10); % trajectory points colour mapped to heat rate
end
plot(longitude_Plot(1)*180/pi,-latitude_geo(1)*180/pi,'g.','MarkerSize',30); % entry interface point
text(0.98*longitude_Plot(1)*180/pi,-0.98*latitude_geo(1)*180/pi,'Entry Interface (120 km)');
plot(landingLong,-landingLat_geo,'r.','MarkerSize',20); % landing site point
text(0.98*landingLong,-0.98*landingLat_geo,'Landing Site');
axis equal
xticks([-180:30:180]);
yticks([-90:15:90]);
hax = gca;
hax.YTickLabel = flipud(hax.YTickLabel);
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longtiude (°)');
ylabel('Latitude (°)');
grid on

%% Stop Clock
fprintf('\n');
toc

function dzdt = dz(t,z,G,M,BC,R0,Cl,Cd,wEarth,f)
% z1 = x, z2 = y, z3 = z,
% z4 = vx, z5 = vy, z6 = vz
[long,lat,r] = cart2sph(z(1), z(2), z(3));
v_Air = [z(4), z(5), z(6)] - wEarth*norm(r)*[-sin(long)*cos(lat), cos(long)*cos(lat), 0]; % air speed
altitude = r - R0^3*(1-f)/sqrt(R0^4*(1-f)^2*(cos(long)^2+sin(long)^2)*cos(lat)^2+R0^4*sin(lat)^2);

aG = G*M/(norm(r)^2); % gravity
aD = EarthAtmosRho(altitude,R0)/2/BC*norm(v_Air)^2; % drag
aL = aD/Cd*Cl; % lift

h = cross([z(1), z(2), z(3)],v_Air); % angular momentum vector
liftVect = cross(v_Air,h); % lift vector (up)

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

function [res] = horizontalViewDist(x,R0,h)
res = x^2+(R0^2-x^2)-(R0+h)*sqrt(R0^2-x^2);
end
