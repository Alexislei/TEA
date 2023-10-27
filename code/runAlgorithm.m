clear all
clc

% Configure any of these parameters to match needs of given problems

popsize=param(:,1); % population size P
Gen=param(:,2); % generation size
R=param(:,3); % the percentage of replication of well-performed chromosomes  5%
pc=param(:,4); % crossover rate
pml=param(:,5); % mutation rate local
pmg=param(:,6); % mutation rate global


load data_kazeroon.mat
Ge=2000;
params=[200 Ge 0.04 0.6 0.08 0.01];
[cost_best,pop_per_gen]=main_func_n1_1(params,F,C);
save(['n1_1_',num2str(double(a)),'.mat'],'cost_best','pop_per_gen');

