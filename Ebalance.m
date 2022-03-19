clear all
close all
clc

H_ads = 18.6e3; % mean heat of adsorption J/mol 
cp_air = 29.13; % heat capacity air J/mol K
cp_cms = (8.598+9.251)/2; % heat capacity cms J/mol K

e = 0.36; % voidage
D = 0.65; % diameter m
H = 2.1; % height m
V_cms = 0.25*pi*D^2*H*(1-e);
rho = 1000; % density kg/m3
m_cms = V_cms*rho; %kg
MW = 12e-3; % kg/mol
t_cycle = 22;

Nin = 167.02; % inlet mol flowrate mol/s
Tin = 20+273; 

t = linspace(0,t_cycle, 100);

[t, T] = ode45(@f,t,[Tin], [], cp_cms, m_cms, cp_air, Nin, Tin, H_ads, MW);

Nout = 0.9728.*t+17.83; % mol/s
Qcool = (Nin-Nout).*(cp_air*Tin+H_ads)/1000; % kW

tiledlayout(1,2)
nexttile
plot(t, T)
xlabel('Time (s)')
ylabel('Temperature (K)')
xlim([0 t_cycle])
title('a)')
nexttile
plot(t, Qcool)
xlabel('Time (s)')
ylabel('Cooling Duty (kW)')
xlim([0 t_cycle])
title('b)')

function dTdt=f(t, T, cp_cms, m_cms, cp_air, Nin, Tin, H_ads, MW)

    Nout = 0.9723*t+18.1493; % mol/s
    N = 137.1314*t+45.5463; % mol
    dNdt = 137.187;

    tau = cp_cms*m_cms/MW+cp_air*N;

    dTdt= (1/tau)*(cp_air*Nin*Tin + Nin*H_ads - Nout*H_ads - cp_air*T*(Nout+dNdt));
end
