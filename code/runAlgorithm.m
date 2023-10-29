clear all
clc

% Configure any of these parameters to match the needs of given problem
params.gen=gen; % number of generations
params.popsize=popsize; % population size
params.m=m; params.n=n; % individual size
params.pc=pc; % crossover rate
params.pml=pml; % mutation rate local
params.pmg=pmg; % mutation rate global

data=load("data_kazeroon.mat");
[fit_best,pop_best]=func_TEA(params,data); % require parameters and any problem-dependent data
save(['result_TEA','.mat'],'fit_best','pop_best','-v7.3');

