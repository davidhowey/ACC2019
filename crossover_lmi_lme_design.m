% Crossover system, Adaptive Observer
% Solution LMI - LME Problem
% Mosek and Yalmip needed

clear all
close all

ops = sdpsettings('verbose',1,'savesolveroutput',1,'showprogress',1,...
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


%Space State Matrices 
A=[ 0 0; dot_V/(epsil*V_cell), -dot_V/(epsil*V_cell)];
B=-(1/c_0)*[1/(Far*V_res); 1/(epsil*Far*V_cell)];
E=-(1/c_0)*[1/V_res; 1/(epsil*V_cell)];
C=[0 1];

%Observability
%[C; C*A]
Ob1=obsv(A,C); 

cond(Ob1)
%Number of unobservable states
unob = length(A)-rank(Ob1)


%Adaptive Observer Design
[n,n]=size(A);
[n,m]=size(E);
[q,n]=size(C);

P = sdpvar(n,n,'symmetric');
W = sdpvar(n,n);
Z= sdpvar(n,q);
F= sdpvar(m,q);

CT=[...
    P>=0,...
    W>=eye(2),...
    A'*P+P*A-C'*Z'-Z*C+[1 0; 0 100]*W<=0,...
    E'*P-F*C==0,...
    ];

sol=optimize(CT,[],ops);
P_v=value(P);
W_v=value(W);
Z_v=value(Z);
F_v=value(F);

eig(P_v)
eig(W_v)
eig(A'*P_v+P_v*A-C'*Z_v'-Z_v*C+[1 0; 0 100]*W_v)
L_v=inv(P_v)*Z_v
F_v

save('obs_gains', 'L_v', 'F_v');





 
