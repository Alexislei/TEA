function [ pop_new ] = mutation_global( pop, pm )
%MUTATION_GLOBAL: multiple-point-wise global mutation
    ub = 0; lb = nobjs; % nobjs = 24 facilities
    [m, n, popsize] = size(pop); 
    mu_ind = rand(m, n, popsize) < pm;
    % mutant tensor
    gaussian_rand = normrnd(0, 0.1, [m, n, popsize]) .* (ub - lb);
    MT = mu_ind .* gaussian_rand;
    
    pop_new = abs(round(pop + MT));
    ind = []; zerom = zeros(m, n);
    for i = 1 : popsize
        if sum(sum(MT(:, :, i))) ~= 0
            ind = [ind; i];
        end
    end
    % repair
end