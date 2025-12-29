%% AME 30335 - Heat Transfer
%  Electronic Cooling System Design Project 
%  Part 2: Cold Plate Design, q''_required = 100 W/cm^2
%  Anya Lindholm & Joe Brady
%  Due Friday, April 26

%% Start
clear; clc; close all;
 
%% Set Parameters
N = 400;
vel= 1;%linspace(0, 2, N); %[m/s]
rho = 997;
D = linspace(0.001, 0.009, N);%[0.004 0.006 0.01 ];%[0.008];
L = 0.20;%[02020.32 0.20 0.12]; %[m]
%D = linspace(L/50, L/5, 200);
%D = 0.006;
Pr = 6.62;
Pthic = 0.0005;
kWater = 0.606; %[W/mK]
mu = 959e-6; %[N*s/m^2]
kCopper = 401; %[W/mK]
Ts = 90; %[Celsius]
Tmi = 20; %[Celsius]
cp = 4181;

%% Calculate flux and pressure drop
for d = 1:length(D)
    for i = 1:length(vel)
        %for l = 1:length(L)
        l = 1;
        %find the reynolds number
        A = pi.*(D(d)/2).^2; %[m^2]
        mdot = rho.*vel(i).*A; %[m]
        Re = (4.*mdot)./(D(d).*pi.*mu);
        xh = 0.05*Re*D(d);
        xth = Pr*xh;
        % see if flow is fully developed
        if L(l)/D(d) >= 10 || xth < L(l)
            % Determine Nu based on Reynolds Number
            if Re < 2300
                Nu = 3.66;
            end
            if Re > 10000
                Nu = 0.023*(Re^(0.8))*Pr^(0.4);
                cond = 'turbulent';
            end
            if Re >= 2300 && Re <= 10000
                f = (0.79*log(Re)-1.64).^(-2);
                Nu = ((f/8).*(Re-1000)*Pr) ./ (1 + 12.7*((f/8).^0.5)*(Pr^(2/3)-1)); 
                cond = 'in between';
            end
        else 
            Nu = 3.66 + 0.0668*(D(d)/L(l))*Re*Pr/(1+0.04*((D(d)/L(l)*Re*Pr)^(2/3)));
            cond = 'developing';
        end

        
        
        %find h and from there 
        h = Nu.*kWater./D(d);
        Rtot = log(1+Pthic)/(2*pi*L(l)*kCopper)+1/(h*pi*D(d)*L(l));
        Tmo = Ts - exp(-1./(mdot.*cp.*Rtot))*(Ts-Tmi);
        deltaT = ((Ts-Tmo) - (Ts-Tmi)) ./ log((Ts-Tmo) ./ (Ts-Tmi));
        q(d) = deltaT/(Rtot)/16;
        %flux(d,i) = q(d,i)/(0.4^2);
        f = (0.79*log(Re)-1.64).^(-2);
        deltaP(d) = f.*rho.*vel(i).^2 ./ (2.*D(d))*L(l);
        Ppump(d) = deltaP(d)*mdot/rho;
    end
    end
%end

%% plot results
figure(1)
hold on
plot(100*D, q(1,:))
% plot(vel, q(2,:))
% plot(vel, q(3,:))
% plot(vel, q(4,:))
xlabel('diameter, D (cm)', 'interpreter', 'latex', 'fontsize',16)
ylabel('Cooling Heat Flux ($\frac{W}{cm^2}$)', 'interpreter', 'latex', 'fontsize',16)
%legend('D = 0.4cm, 8 passes', 'D = 0.6cm, 5 passes',  'D = 1cm, 3 passes', 'fontsize', 11)
grid on
hold off

figure(2)
hold on
    plot(100*D, deltaP(1,:))
%     plot(vel, deltaP(2,:))
%     plot(vel, deltaP(3,:))
%     plot(vel, deltaP(4,:))
xlabel('diameter, D (cm)', 'interpreter', 'latex', 'fontsize',16)
ylabel('pressure drop, $\Delta p$ (Pa)', 'interpreter', 'latex', 'fontsize',16)
%legend('D = 0.4cm, 8 passes', 'D = 0.6cm, 5 passes',  'D = 1cm, 3 passes', 'fontsize', 11)
grid on
hold off

 figure(3)
    plot(D, Ppump(1,:))
%     plot(vel, deltaP(2,:))
%     plot(vel, deltaP(3,:))
%     plot(vel, deltaP(4,:))
xlabel('fluid velocity, u (m/s)', 'interpreter', 'latex')
ylabel('pressure drop, $\Delta p$ (Pa)', 'interpreter', 'latex')
%legend('D = 0.4cm, 8 passes', 'D = 0.6cm, 5 passes',  'D = 1cm, 3 passes')
grid on
hold off
