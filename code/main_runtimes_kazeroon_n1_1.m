clear all;clc;
runtimes=30;
cost_be=zeros(runtimes,1);

for a=1:runtimes
    load data_kazeroon.mat
    Ge=2000;
    params=[200 Ge 0.04 0.6 0.08 0.01];
    [cost_best,pop_per_gen]=main_func_n1_1(params,F,C);
    save(['n1_1_',num2str(double(a)),'.mat'],'cost_best','pop_per_gen');
end