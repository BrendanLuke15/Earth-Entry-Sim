# Earth Entry Simulator
**This is a relatively simple low-fidelity sim for blunt capsules (re)entering Earth's atmosphere.**

Some features:
- Simple atmospheric model from [Braeunig's Rocket & Space Technology](http://www.braeunig.us/space/atmmodel.htm) (US Standard Atmos)
- 2D & 3D versions
- fixed aerodynamic coefficients
- 1st order entry heating and heatshield design (PICA)

-----

**Vehicle Design:**

Configure the mass, geometry, and aerodynamics of the vehicle (point mass). SpaceX Crew Dragon example from my work [here](https://space.stackexchange.com/a/55685/40257):

```
% input values: DESIGN
m = 21200/2.205; % mass of entry vehicle (kg)
Cd = 1.23; % drag coefficient
Cl = 0.27*Cd; % lift coefficient                    -> modulated in flight (roll angle), used as constant "fudge factor" in this sim
d = 4; % entry vehicle diameter (m)
S = (d/2)^2*pi(); % ref. area (m^2)
BC = m/(S*Cd); % ballistic coefficient (kg/m^2)
Rn = 2; % nose radius (m)                           -> used as a "fudge factor" in this sim for heating calcs
emis = 0.9; % thermal IR emmisivity (0.9 for PICA)
```

**Entry Conditions Setup:**

Pseudo orbital elements: apoapsis height, peripsis height, inclination, (Earth) longitude of ascending node, argument of peripsis. Again SpaceX Crew Dragon example [here](https://space.stackexchange.com/a/55685/40257) & [here](https://space.stackexchange.com/a/58332/40257):

```
Conditions = [22,418]; % periapsis and apoapsis altitudes (km)

% Initial Orbit Angles
inc = 51.6; % inclination (°)                            -> ISS inclination
RAAN = 127; % (Earth) longitude of ascending node (°)    -> RAAN & w chosen to facilitate landing near Florida from the northwest
w = 165; % argument of periapsis (°)
```

**Entry Heating & PICA Heat Shield Design:**

Estensive leverage of methods from [my work here](https://space.stackexchange.com/a/55725/40257). Returns a first order estimate of required stagnation point 1D heatshield thickness. Specific to the PICA material:
- Calculate Sutton & Graves convective & Martin radiative heating rates
- assume radiative equilibrium to find thermal environment around vehicle stagnation point
- use PICA ablation data and thermal properties to:
  - estimate & integrate ablation rate to find total thickness "lost"
  - 1D unsteady heat flow to determine minimum "insulation thickness" to keep bondline from failing


Heating & trajectory data can also be used to plan observations of 
