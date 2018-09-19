% Cross over system
% Data
% 10-07-2018

clear all
close all

v_eq = 2.2;
eq_range = 0.4;

xover_data_full = readall(tabularTextDatastore(strcat('data.txt')));
xover_data_full.(4) = xover_data_full.(4)*3600;
xover_data = xover_data_full(xover_data_full.(3)==9,:);
xover_data.(4)=xover_data.(4)-xover_data.(4)(1);
xover_data = xover_data(xover_data.(9)>v_eq-eq_range/2,:);
xover_data= xover_data(xover_data.(9)<v_eq+eq_range/2,:);
v_data = xover_data.(9);
t_data = xover_data.(4)-xover_data.(4)(1)+0.001;

plot (t_data, v_data);
save('crossover_data', 't_data', 'v_data');
