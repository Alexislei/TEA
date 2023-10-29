function [ pop ] = init( popsize, m, n)
%INIT: initializes the population
    pop=zeros([m, n, popsize]);
    
    % nobjs = 24; % number of facilities
    
    for i = 1 : popsize
        chrom = pop(:,:,i);
        % initialize each individual based on the given problem
        for k = 1 : nobjs
            x = unidrnd(m); y = unidrnd(n); % random position (x,y)
            while chrom(x,y) ~= 0 % skip occupied position
                x = unidrnd(m);y = unidrnd(n);
            end
            chrom(x,y) = k;
        end
        pop(:, :, i) = chrom;
    end
end