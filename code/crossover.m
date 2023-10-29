function [ pop_new ] = crossover( pop, pc )
% CROSSOVER: T-product crossover
    [m, n, popsize] = size(pop); 
    pop_new = zeros(size(pop));
     
    for k = 1 : 2 : (popsize-1)
%         cid = randperm(popsize, 2);
%         cid=[k:k+1];
        X = pop(:, :, cid);
        
        if rand < pc
            % strating point of selected rows & columns
            sp_r = unidrnd(m); sp_c = unidrnd(n);
            % number of selected rows & columns [1,m]&[1,n]
            L_r = unidrnd(m); L_c = unidrnd(n);
            diag_r = ones(m,1); diag_c = ones(n,1);
              
            for i = 1 : L_r
                if mod(sp_r + i-1, m) == 0
                    diag_r(m,:) = 0;
                else
                    diag_r(mod(sp_r + i - 1, m), :) = 0;
                end
            end
            
            for i = 1 : L_c
                if mod(sp_c + i-1, n) == 0
                    diag_c(n, :) = 0;
                else
                    diag_c(mod(sp_c + i - 1, n), :) = 0;
                end
            end
            
            A(:, :, 1) = diag(diag_r); A(:, :, 2) = diag(1 - diag_r);
            B(:, :, 1) = diag(diag_c); B(:, :, 2) = diag(1 - diag_c);
            
            if rand < 0.5
                X_new = tprod(A, X);
            else
                X_new = tprod(X, B);
            end
            pop_new(:, :, k : k+1) = X_new;
        else
            pop_new(:, :, k : k+1) = X;
        end
    end
    
end