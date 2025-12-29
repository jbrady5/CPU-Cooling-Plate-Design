%% AME 30335 - Heat Transfer
%  Electronic Cooling System Design Project
%  Anya Lindholm & Joe Brady
%  Due Friday, April 26

% Design two heat sinks for two CPUs to achieve 1 W/cm^2

%% Start
clear; clc; close all;

%% 7.6 - Flow Across Banks of Tubes 
% Parameters
% Vary length L from 0.5 cm to 4 cm
Larray = linspace(0.005, 0.04, 50);

for (i=1:50)
% Chip Parameters
T = 90; %[Celcius]
L = 0.04; %[m]
A = L^2; %[m^2]

Tinf = 23; %[Celcius]
Ti = Tinf; %[Celcius]
Vinf = 0.1; %[m/s]

% Pin Parameters
D = 0.006; % Diameter [m]
P = pi*D; % Perimeter [m]
As = (pi/4)*D^2; % Surface area [m^2]

% Array Parameters
N = 50; % Number of tubes
N_L = 10; % Number of tube rows
Dtubes = (0.03 - 4*D)/3; % Distance between tubes in same row
S_T = Dtubes + D;
S_D = S_T;
S_L = sqrt(S_D^2 - (S_T/2)^2);
PL = S_L/D;
PT = S_T/D;

% Calculated Parameters
Theta = T - Tinf;
% Eqns. 7.60 & 7.61
if (2*(S_D-D) >= (S_T-D))
    Vmax = S_T*Vinf/(S_T-D);
end
if (2*(S_D-D) < (S_T-D))
    Vmax = S_T*Vinf/(2*(S_D-D));
end

%% Row Temperature Assumptions
% Assume a 1 degree celcius drop after each row
for (j=1:N_L)
    T(j) = Ti + (j-1);
    TAvg(j) = (Ti+T(j))/2;
end

% Evaluate properties at each location interpolation valid for -23C < Tinf < 27C   
for (j=1:N_L)
     % Evaluate at Tavg  (kg/m^3)
      rho(j) = ((TAvg(j)+273)-250)/(300-250)*(1.1614-1.3947)+1.3947; 
     % Evaluate at Tavg  (N*s/m^2)
      mu(j) = (((TAvg(j)+273)-250)/(300-250)*(184.6-159.6)+159.6)*10^-7;
     % Evaluate at Tavg  (W/m*K)
      k_air(j) = (((TAvg(j)+273)-250)/(300-250)*(26.3-22.3)+22.3)*10^-3;
     % Evaluate at T
      Pr(j) = ((T(j)+273)-250)/(300-250)*(0.707-0.720)+0.720;
end
    
% Values for Cp are all pretty similar
Cp_air = 1.007E-3; %J/kg*K

%% Check Pressure Drop
X = 1; %Table 7.15
f = 1;
deltaP = N_L*X*((rho(1)*Vmax^2)/2)*f; % [Pa]

%% Nusselt Number and Convection Coefficient
for (j=1:N_L)
    Re(j) = rho(j)*Vmax*D/mu(j);
end

% Evaluate constants C1, C2, m using Table 7.5 and 7.6
% Shown valid for 10^2 < Re < 10^3 (treat like single cylinder)

for (j=1:N_L)
    Nu(j) = 0.3+0.62*Re(j)^(1/2)*Pr(j)^(1/3)*(1+(Re(j)/282000)^(5/8))^(4/5)/((1+(0.4/Pr(j))^(2/3))^(1/4));
    h(j) = Nu(j)*k_air(j)/D;
end

%% Fin Analysis
% Aluminum
k_fin = 237; % Aluminum (makes the graph highkey gross)
% k_fin = 1.4 % (makes the graph the correct shape)

%hAvg = mean(h);
m = sqrt(h(j)*P/(k_fin*As));

for (j=1:N_L)
    M(j) = sqrt(h(j)*P*k_fin*As)*Theta;
    Qf(j) = M(j)*(sinh(m*Larray(i))+(h(j)/(m*k_fin))*cosh(m*Larray(i)))./(cosh(m*Larray(i))+(h(j)/(m*k_fin))*sinh(m*Larray(i)));
end

% Note: Qtot dependent on fin array geometry!
Qtotal = Qf(1)*4 + Qf(2)*3 + Qf(3)*4 + Qf(4)*3; % [Watts]
Flux(i) = (Qtotal/A)/(100^2); % [Watts/cm^2]

end

%% Plot
figure(1)
plot(Larray*100,Flux,'b')
xlabel('Fin length (cm)', 'interpreter', 'latex')
ylabel('Cooling Heat Flux ($\frac{W}{cm^2}$)', 'interpreter', 'latex')
