% Crossover system, Adaptive Observer
% Solution LMI - LME Polytopic Problem
% Mosek and Yalmip needed


clear all
close all

ops = sdpsettings('verbose',1,'savesolveroutput',1,'showprogress',0,...
    'solver','mosek-sdp',...
     'mosek.MSK_DPAR_INTPNT_TOL_PFEAS',1e-11,...
     'mosek.MSK_DPAR_INTPNT_TOL_DFEAS',1e-11,...
     'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP',1e-11);

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
A=@(q) [ 0 0; q/(epsil*V_cell), -q/(epsil*V_cell)];
B=-(1/c_0)*[1/(Far*V_res); 1/(epsil*Far*V_cell)];
E=-(1/c_0)*[1/V_res; 1/(epsil*V_cell)];
C=[0 1];

Qdom=[0.9*dot_V, dot_V, 1.1*dot_V];


%Observability
%[C; C*A]
Ob1=obsv(A(Qdom(1)),C); 
Ob2=obsv(A(Qdom(2)),C); 
Ob3=obsv(A(Qdom(3)),C); 

cond(Ob1)
cond(Ob2)
cond(Ob3)
%Number of unobservable states
unob = length(A(Qdom(1)))-rank(Ob1)


%Adaptive Observer Design
[n,n]=size(A);
[n,m]=size(E);
[r,n]=size(C);

P = sdpvar(n,n,'symmetric');
W = sdpvar(n,n);
Z= sdpvar(n,r);
F= sdpvar(m,r);
balpha=sdpvar(1,1);
gamma_z=sdpvar(1,1);
gamma_f=sdpvar(1,1);
I2=eye(2);
beta=1e-4;


CT=[...
    P>=0,...
    W>=[1e-1, 0; 0 1e-6],...
    P*E-(C'*F')==0,...
    [gamma_z*I2  Z; Z' gamma_z]>=0,...
    [gamma_f  F; F' gamma_f]>=0,...
    gamma_f>=0,...
    gamma_z>=0,...
    balpha>=1e-3,...
    ];

for i=1:length(Qdom)
    
       CT=[CT,...
           [-A(Qdom(i))'*P-P*A(Qdom(i))+C'*Z'+Z*C-(beta*I2+W), P; ...
           P, balpha*I2]>=0,...
           ];
           
end

kappa_f=1e-5;
kappa_z=1;

sol=optimize(CT,[balpha+kappa_z*gamma_z+kappa_f*gamma_f],ops);
P_v=value(P);
W_v=value(W);
Z_v=value(Z);
F_v=value(F);
balpha_v=value(balpha);

eig(P_v)
eig(W_v)
An=A(Qdom(2));
eig(An'*P_v+P_v*An-C'*Z_v'-Z_v*C+ (1/balpha_v)*P_v*P_v+(beta*I2+W_v))
L_v=inv(P_v)*Z_v
F_v



save('obs_gains_poly', 'L_v', 'F_v');




 
