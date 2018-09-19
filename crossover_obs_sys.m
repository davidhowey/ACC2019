% Crossover system
% System Simulation and Adaptive Observer 
% LMI-LME (polytopic) design

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

%Functions
I=@(t) 0*t;
dot_Nx=@(z) k_mt*c_0*(z);                  %[L/s] (slope = -k_mt)


%Space State Matrices 
A=[ 0 0; dot_V/(epsil*V_cell), -dot_V/(epsil*V_cell)];
B=-(1/c_0)*[1/(Far*V_res); 1/(epsil*Far*V_cell)];
E=-(1/c_0)*[1/V_res; 1/(epsil*V_cell)];
C=[0 1];

%Observer gains
%'L_v', 'F_v'
%load obs_gains.mat

%polytopic
load obs_gains_poly.mat

%Simulation
%---------------------------------------------------------------


%---------------------------------------------------------------
%Cross over system:
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
Vout=E0_cell+(R*T*2/Far)*log(SOC_cell./(1-SOC_cell));

%---------------------------------------------------------------
 

%---------------------------------------------------------------
%observer
ltime=length(tout)

%observer initial condition
xe0=[0.85; 0.8];

%parameters
p=7;
th0=0*randn(p,1);
Lambda=0.5*(1e-5)*eye(p)*abs(1/F_v);
sigma=1e-1;

xini=[xe0;th0];

%---------------------------------------------------
%---------------------------------------------------
%basis
Cx=linspace(0.05,0.95,p);
di=[Cx(2)-Cx(1)];
si=1./(0.6*di).^2;
Sx=ones(1,p).*si(1);

%----------------------------------------
%z: 
%row: state variable, 
%columns: time
mon=@(z,i) exp(-0.5*diag((z-Cx(i))'*[(z-Cx(i)).*Sx(i)]));

vs='monv=@(z)[mon(z,1)';
for j=2:p
    vs=[vs,',', 'mon(z,',int2str(j),')'];
end
vs=[vs,']'];

%creating 
%monv=@(z)[mon(z,1),...,mon(z,p)]
eval(vs);
%-----------------------------------------
%result
%rows: 
%columns 

phi=@(z) monv(z)./(sum(monv(z)'))';

%ploting results
figure(2);
zv=(0:0.005:1);
plot(zv,phi(zv)); drawnow
%---------------------------------------------------
%---------------------------------------------------


%System: Process and Adpative Observer
obs_cross_sys=@(t,x,y) [A*x(1:2)+ E*(phi(y)*x(3:end))+ B*I(t)+ L_v*(y-x(2));...
                        Lambda*([phi(y)]'.*(F_v*(y-x(2)))-0.5*sigma*x(3:end)*norm((y-x(2))))];

                  
%Step simulation

SOCe(1)=xini(1);
SOCe_cell(1)=xini(2);
the_cell(:,1)=xini(3:end);

y=SOC_cell(1);
dot_Nxe(1)=phi(y)*the_cell(:,1);

for k=1:ltime-1,

k
tspank=linspace(tout(k),tout(k+1),3);

y=SOC_cell(k+1);

[toute,xsole] = ode23s(@(t,x) obs_cross_sys(t,x,y),tspank,xini);

SOCe(k+1)=xsole(end,1);
SOCe_cell(k+1)=xsole(end,2);
the_cell(:,k+1)=xsole(end,3:end)';
dot_Nxe(k+1)=phi(y)*the_cell(:,k+1);

xini(1,1)=SOCe(k+1);
xini(2,1)=SOCe_cell(k+1);
xini(3:end,1)=the_cell(:,k+1);

Vout_e=E0_cell+(R*T*2/Far)*log(SOCe_cell(1:k+1)./(1-SOCe_cell(1:k+1)));


%Figures
figure(1)

subplot(231);
plot(tout(1:k+1)/3600,SOC(1:k+1),'g',tout(1:k+1)/3600,SOCe,'r--','LineWidth',2); 
title('SOC');
xlabel('Time[hrs]');drawnow

subplot(232);
plot(tout(1:k+1)/3600,the_cell(:,1:k+1));
title('Parameters');
xlabel('Time[hrs]');drawnow

subplot(233);
plot(tout(1:k+1)/3600,SOC_cell(1:k+1),'g',tout(1:k+1)/3600,SOCe_cell,'r--','LineWidth',2);  
title('SOC_{cell}');
xlabel('Time[hrs]');drawnow


subplot(223);
plot(SOC_cell(1:k+1), dot_Nx(SOC_cell(1:k+1)),'g',SOC_cell(1:k+1),dot_Nxe(1:k+1),'r--','LineWidth',2);
title('Cross-over Rate');
xlabel('SOC_{cell}'); drawnow

subplot(224);
plot(SOC_cell(1:k+1), Vout(1:k+1),'g',SOC_cell(1:k+1),Vout_e(1:k+1),'r--','LineWidth',2);
title('Voltage: V_{out} wrt. SOC_{cell}');
xlabel('SOC_{cell}'); drawnow

end

%save crossover_obs_sys.mat
save crossover_obs_sys_poly.mat


 
