% Brendan Luke
% April 24, 2021
clear, clc
format compact
close all
tic
%% Input Parameters
% input values: CONSTANT
R0 = 6378137; % radius of Earth (m)
G = 6.6743*10^-11; % gravitational constant (N*m^2/kg^2)
M = 5.97237*10^24; % mass of Earth (kg)

% input values: DESIGN
m = 9615; % mass of entry vehicle (kg)
Cd = 1.23; % drag coefficient
Cl = 0.13*Cd;%0.07; % lift coefficient
d = 4; % entry vehicle diameter (m)
S = (d/2)^2*pi(); % ref. area (m^2)
BC = m/(S*Cd); % ballistic coefficient (kg/m^2)
Rn = 4; % nose radius (m)
emis = 0.9; % thermal IR emmisivity (0.9 for PICA)
type = 'Orbita';%'Orbital'; % specify type of entry: 'Orbital' or 'IC' (initial conditions)
if strcmp(type,'Orbital')
    Conditions = [22,419]; % periapsis and apoapsis altitudes (km)
else % initial conditions case
%     Conditions = [-5.2,11.4]; % EFPA (°) and Entry Velocity (km/s)
%     Conditions = [-7.8,11.4]; % EFPA (°) and Entry Velocity (km/s)
    Conditions = [-6.45,11.4]; % EFPA (°) and Entry Velocity (km/s)
end

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
if strcmp(type,'Orbital')
    a = 1000*mean(Conditions)+R0; % orbit semi major axis (m)
    V_Entry = sqrt(G*M*(2/(120000+R0)-1/a)); % entry velocity (m/s)
    angMomentum = (Conditions(1)*1000+R0)*sqrt(G*M*(2/(Conditions(1)*1000+R0)-1/a)); % angular momentum m^2/s
    p = angMomentum^2/(G*M); % semi-latus rectum (m)
    e = (1000*diff(Conditions))/(2*R0+1000*sum(Conditions)); % orbital eccentricity
    EFPA = -acosd(angMomentum/V_Entry/(120000+R0)); % entry flight path angle (°)
    fprintf('The entry velocity is %.0f m/s and the EFPA is %.2f°.\n',V_Entry,EFPA)
    
    % create data for plot
    phi = 90-acosd(1/e*(p/(R0+1000*120)-1)); % phase angle to make entry point occur at theta = 90°
    if -sind(phi) > 0
        adj = phi;
    elseif cosd(phi) < 0
        adj = phi-90;
    else
        adj = phi;
    end
    
    for i = 1:361
        r_plot(i) = p/(1+e*cosd(i-phi));
        x_plot(i) = r_plot(i)*cosd(i);
        y_plot(i) = r_plot(i)*sind(i);
    end
else % initial conditions case
    angMomentum  = (R0+120000)*1000*Conditions(2)*cosd(Conditions(1)); % angular momentum (m^2/s)
    p = angMomentum ^2/G/M; % Semi-latus rectum (m)    
    a = 1/(2/(R0+120000)-(1000*Conditions(2))^2/G/M); % orbit semi-major axis (m)
    Vinf = sqrt(-G*M/a); % hyperbolic excess velocity (m/s)
    b = angMomentum/Vinf; % orbit semi-minor axis (m)
    e = sqrt(1+b^2/a^2); % orbital eccentricity
    rp = a*(1-e); % periapsis distance (m)
    fprintf('The entry velocity is %.0f m/s and the EFPA is %.2f°.\n',Conditions(2)*1000,Conditions(1))
    fprintf('The close approach altitude is %.1f km.\n',(rp-R0)/1000)
    
    % create data for plot
    phi = 90-acosd(1/e*(p/(R0+1000*120)-1)); % phase angle to make entry point occur at theta = 90°
    if -sind(phi) > 0
        adj = phi;
    elseif cosd(phi) < 0
        adj = phi-90;
    else
        adj = phi;
    end
    
    for i = 1:361
        r_plot(i) = p/(1+e*cosd(i-phi));
        x_plot(i) = r_plot(i)*cosd(i);
        y_plot(i) = r_plot(i)*sind(i);
    end
end

% Constant Things
for i = 1:361
    % fill in earth
    x_fill(i) = R0*cosd(i);
    y_fill(i) = R0*sind(i);
    % atmosphere line
    xAtmos(i) = (120+R0/1000)*cosd(i);
    yAtmos(i) = (120+R0/1000)*sind(i);
end

% make plot
figure('name','Entry Geometry');
set(gcf,'WindowState','maximized');
fill(x_fill/1000,y_fill/1000,[0, 0.55, 0.8],'FaceAlpha',0.4);
hold on
if strcmp(type,'Orbital')
    plot(x_plot/1000,y_plot/1000,'LineWidth',2,'Color',[1 0 0]);
    plot((R0/1000+Conditions(2))*cosd(180+adj),(R0/1000+Conditions(2))*sind(180+adj),'k.','MarkerSize',20); % apoapsis point
    plot((R0/1000+Conditions(1))*cosd(adj),(R0/1000+Conditions(1))*sind(adj),'k.','MarkerSize',20); % periapsis point
    text(1.08*(R0/1000+Conditions(2))*cosd(180+adj),1.08*(R0/1000+Conditions(2))*sind(180+adj),['Apoapsis ',num2str(Conditions(2)),' km'],'HorizontalAlignment','Right');
    text(1.08*(R0/1000+Conditions(1))*cosd(adj),1.08*(R0/1000+Conditions(1))*sind(adj),['Periapsis ',num2str(Conditions(1)),' km']);
    plot(xAtmos,yAtmos,'k-');
    plot(0,R0/1000+120,'g.','MarkerSize',30); % entry interface point
    text(0,1.08*(R0/1000+120),'Entry Interface');
else % initial conditions case
    plot(x_plot/1000,y_plot/1000,'LineWidth',2,'Color',[1 0 0]);
    plot(xAtmos,yAtmos,'k-');
    plot(0,R0/1000+120,'g.','MarkerSize',30); % entry interface point
    text(0,1.08*(R0/1000+120),'Entry Interface');
end

grid on
axis equal
xlabel('X (km)');
ylabel('Y (km)');

if strcmp(type,'Orbital')
    axis([-1.5*R0/1000 1.5*R0/1000 -1.5*R0/1000 1.5*R0/1000]);
else % initial conditions case
    axis([-1.5*R0/1000 1.5*R0/1000 -1.5*R0/1000 1.5*R0/1000]);
end

%% Trajectory Simulation
y0 = 120000+R0; % initial height
x0 = 0; % initial horizontal distance

if strcmp(type,'Orbital')
    z0 = [x0; y0; V_Entry; EFPA/180*pi]; 
else
    z0 = [x0; y0; 1000*Conditions(2); Conditions(1)/180*pi]; 
end

ReEntryEventsFcn2 = @(t,z) ReEntryEventsFcn(t,z,R0);
options = odeset('MaxStep',100,'Stats','off','RelTol',1e-9,'Events',ReEntryEventsFcn2);
[t,z] = ode45(@(t,z) dz(t,z,G,M,BC,R0,Cl,Cd),[0 100*86400],z0,options);

%% Process Data
altitude = z(:,2)-R0; % meters
for i = 1:length(t)
    temp = EarthAtmos(z(i,2)-R0,R0);
    temp2 = dz(t(i),z(i,:),G,M,BC,R0,Cl,Cd);
    aG(i) = G*M/(z(i,2).^2);
    aD(i) = EarthAtmosRho(z(i,2)-R0,R0)/2/BC*z(i,3)^2;
    Q(i) = 1/2*temp(3)*z(i,3)^2;
    Pt(i) = Q(i) + temp(1);
    Mach(i) = z(i,3)/temp(5);
    qConv(i) = 1.7415*10^-4*sqrt(temp(3)/Rn)*z(i,3)^3/100/100;
    qRad(i) = 8*2.787*10^-67*Rn^0.2*temp(3)^1.05*z(i,3)^19/100/100;
    qTot(i) = qConv(i) + qRad(i);
    vertAccel(i) = temp2(3)*sin(z(i,4)) + z(i,3)*temp2(4)*cos(z(i,4));
    horiAccel(i) = temp2(3)*cos(z(i,4)) - z(i,3)*temp2(4)*sin(z(i,4));
    GForce(i) = aD(i)/9.81;
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
    Tw(i) = ((100^2*qTot(i))/(5.67*10^-8*emis) + temp(2)^4)^0.25;
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

%% Plot Data
figure('name','TPS Performance Map');
set(gcf,'WindowState','maximized');
plot(Q/101325,qTot,'LineWidth',2);
hold on
plot(PICA_Envelope(1:3,1),PICA_Envelope(1:3,2),'k--','HandleVisibility','off','LineWidth',1.5);
plot(PICA_Envelope(4:6,1),PICA_Envelope(4:6,2),'k--','HandleVisibility','off','LineWidth',1.5);
text(0.1,1400,'PICA Envelope');
text(0.3,1650,'PICA Envelope Expansion');
% axis([0 1.5 0 2000]);
xlabel('Stag. Pressure (atm)');
ylabel('Stag. Heat Flux (W/cm^{2})');
grid on

figure('name','G-Force');
set(gcf,'WindowState','maximized');
plot(t,GForce,'LineWidth',2);
grid on
xlabel('Time (s)');
ylabel("G-Froce (g's)");

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

figure('name','Vertical Acceleration');
set(gcf,'WindowState','maximized');
plot(t,vertAccel,'LineWidth',2);
xlabel('Time (s)');
ylabel('Vertical Acceleration (m/s^2)');
grid on

figure('name','Loading Profile');
set(gcf,'WindowState','maximized');
plot(GForce,altitude/1000,'LineWidth',2);
xlabel('G-Force (g)');
ylabel('Altitude (km)');
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

figure('name','Altitude VS Velocity');
set(gcf,'WindowState','maximized');
plot(z(:,3)/1000,altitude/1000,'LineWidth',2);
xlabel('Velocity (km/s)');
ylabel('Altitude (km)');
grid on

figure('name','Simulated Trajectory');
set(gcf,'WindowState','maximized');
plot(xPos/1000,yPos/1000,'LineWidth',2,'Color',[1 0 0]);
hold on
fill(x_fill/1000,y_fill/1000,[0, 0.55, 0.8],'FaceAlpha',0.4);
grid on
axis equal
xlabel('X (km)');
ylabel('Y (km)');
axis([min(xPos)/1000-300 max(xPos)/1000+300 min(yPos)/1000-300 max(yPos)/1000+300]);
% axis([min(xPos)/1000*1.1 max(xPos)/1000*1.1 min(yPos)/1000*1.1 max(yPos)/1000*1.1]);

%% Stop Clock
fprintf('\n');
toc

function dzdt = dz(t,z,G,M,BC,R0,Cl,Cd)
aG = G*M/(z(2)^2);
aD = EarthAtmosRho(z(2)-R0,R0)/2/BC*z(3)^2;
% aD = 0;
% z1 = s, z2 = R, (path dist. and radius)
% z3 = V, z4 = FPA
dzdt(1) = z(3)*cos(z(4));
dzdt(2) = z(3)*sin(z(4));
dzdt(3) = -aD - aG*sin(z(4));
dzdt(4) = aD/z(3)*(Cl/Cd) - (aG - z(3)^2/z(2))/z(3)*cos(z(4));
dzdt = dzdt';
end

% Stop at Skip-out
% function [position,isterminal,direction] = ReEntryEventsFcn(t,z,R0)
% position = [z(2)-R0, z(2)-R0-120000]; % The value that we want to be zero, altitude or skip out case
% isterminal = [1, 1];  % Halt integration 
% direction = [0, 1];   % The zero can be approached from either direction for altitude, increasing for skipout
% end

% Continue at Skip-out
function [position,isterminal,direction] = ReEntryEventsFcn(t,z,R0)
position = [z(2)-R0]; % The value that we want to be zero, altitude or skip out case
isterminal = [1];  % Halt integration 
direction = [0];   % The zero can be approached from either direction for altitude, increasing for skipout
end
