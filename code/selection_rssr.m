function [ pop_new ] = selection_rssr( pop, fitness )
%SELECT: select #popsize chromosomes from population by remainder stochastic sampling with replacement
    
    [m, n, popsize] = size(pop);
    
%     num_expectation = zeros(popsize,1);
    num_expectation = popsize .* (fitness ./ sum(fitness));

%     sum(num_expectation_int)
    num_expectation_int = floor(num_expectation);
    
    pop_new = zeros(m, n, popsize); 
    
    s = 0; % count
    for i = 1 : popsize
        if num_expectation_int(i) ~= 0
            a = [];
            for j = 1 : num_expectation_int(i)
                if j == 1
                    a(:, :, 1) = pop(:, :, i);
                else
                    a(:, :, end+1) = pop(:, :, i);
                end
            end
            pop_new(:, :, s+1:s+num_expectation_int(i)) = a;
            s = s + num_expectation_int(i);
        end
    end
    
   
    % select the remaining undetermined #n-sum(num_expectation_int) individuals for the next generation
    fitness_new = num_expectation - num_expectation_int; 
%     fitness_new=fitness-sum(fitness).*num_expectation_int./popsize;
%     for i=1:popsize
%         fitness_new(i)=fitness(i) - num_expectation_int(i)*sum(fitness)/popsize;%num_expectation=popsize.*(fitness./sum(fitness))
%     end
    
    num_remain = popsize - s;
    fitness_table = cumsum(fitness_new ./ sum(fitness_new));% accumulate fitness_table = zeros(popsize,1);
    pop_remain = [];
    rs = sort(rand(num_remain, 1));
    fiti = 1;
    newi = 1;
    while newi <= num_remain
        if rs(newi) < fitness_table(fiti)
            pop_remain(:, :, newi) = pop(:, :, fiti);
            newi = newi + 1;
        else
            fiti = fiti + 1;
        end
    end
    pop_new(:, :, s+1:end) = pop_remain;
end
