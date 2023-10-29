
function [f_best, p_best] = func_TEA(params, data)

    pop = init(params.popsize, params.m, params.n); % initialization
    fit = evaluate(pop, data);
    
    for t = 1 : params.gen
%         [pop, fit] = re_gen(pop,fit, params.R); % the percentage of replication of well-performed chromosomes  5%
        pop = selection_rssr(pop, fit);
        pop = crossover(pop, params.pc); %crossover 
        pop = mutation_local(pop, params.pml); %mutate
        pop = mutation_global(pop, params.pmg); 
        fit = evaluate(pop,data);
    end
    [f_best, best_idx] = max(fit);
    p_best = pop(:, :, best_idx);
end

% function [ pop_new, fit_new] = re_gen(pop, fitness, R)
% % RE-GEN: remove the worst floor(popsize*R) chromosomes
% % duplicate and insert the best chromosomes into the current population
% 
%     popsize = size(pop, 3);
%     num = floor(popsize * R);
%     pop_new = pop;
%     fit_new = fitness;
%     
%     [~, ind_worst] = sort(fitness, 'ascend'); % worst chromosomes
%     [~, ind_best] = sort(fitness, 'descend'); % best chromosomes
%     
%     pop_new(:, :, ind_worst(1:num, :)) = pop(:, :, ind_best(1:num, :));
%     fit_new(ind_worst(1:num, :), :)=fitness(ind_best(1:num, :), :);
% 
% end
