%Crossover system
%Simulation

clear all
close all

%Parameters
R = 8.314472;       %[J K^-1 mol^-1]
T = 22 + 273;       %[K]
Far = 96485;        %[C/mol]
 
V_res=17.6e-3;      %[L]
c_0=0.1;            %[mol/L]
dot_V= 9e-3/60;     %[L/s]
V_cell = 1.6408*1.3408*(.125*2.54)/(10^3); % volume of one half of reactor chamber in L
epsil=0.87;         %[-]
k_mt=3.3685e-6;     %[L/s] (slope = -k_mt)
E0_cell=2.2;        %[V]  (equilibrium voltage)


%Space State Matrices 
A=[ 0 0; dot_V/(epsil*V_cell), -dot_V/(epsil*V_cell)];
B=-(1/c_0)*[1/(Far*V_res); 1/(epsil*Far*V_cell)];
E=-(1/c_0)*[1/V_res; 1/(epsil*V_cell)];
C=[0 1];

%Functions
I=@(t) 0*t;
dot_Nx=@(z) k_mt*c_0*(z);                  %[L/s] (slope = -k_mt)

%Simulation
%---------------------------------------------------------------


%System:
cross_sys=@(t,x) [A*x(1:2,:)+ E*dot_Nx(x(2,:)) + B*I(t)];


%initial condition
x0=[1;1];  %[SOC, SOC_cell]
t0=0;
te=4.5*3600;                 %4.5 hours
tspan=linspace(t0,te,1000);

%simulation
[tout,xsol] = ode45(@(t,x) cross_sys(t,x),tspan,x0);

%outputs
SOC=xsol(:,1);
SOC_cell=xsol(:,2);

%a=0.9826;
%b=0.03086;
%Vout=E0_cell+(R*T*2/Far)*log((SOC_cell-b)./(a-(SOC_cell-b)));
Vout=E0_cell+(R*T*2/Far)*log(SOC_cell./(1-SOC_cell));


%Figures
figure(1)

subplot(221);
plot(tout/3600,SOC,'b')
title('SOC');
xlabel('Time[hrs]');

subplot(222);
plot(tout/3600,SOC_cell,'r');
title('SOC_{cell}');
xlabel('Time[hrs]');

subplot(223);
plot(tout/3600,Vout,'b');
title('V_{out}: Output Voltage');
xlabel('Time[hrs]');

subplot(224);
plot(SOC_cell,Vout,'r');
title('V_{out} wrt. SOC_{cell}');
xlabel('SOC_{cell}');
ylabel('V_{out}');
 
