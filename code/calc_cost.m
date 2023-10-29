
function [ fit ] = evaluate(pop, data)
% Evaluation function: evaluate the fitness of each individual
    popsize = size(pop, 3);
    fit = zeros(popsize, 1); % fitness vector
    for i = 1 : popsize 
        cost = calc_fit(pop(:, :, I), M); % define the fitness function based on the target problem
        fit(i,:) = cost;
    end
end

