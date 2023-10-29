function [ pop_new ] = mutation_local( pop, pm )
% MUTATION_LOCAL: local t-product mutation
    popsize = size(pop, 3); 
    pop_new = zeros(size(pop));
        
    for i = 1 : popsize 
%         ind = unidrnd(popsize);
        ind = i;
        if rand < pm
            off = abs(round(mutate_convlike(pop(:, :, ind))));
            pop_new(:,:,i)=off;
        else
            pop_new(:,:,i)=pop(:,:,ind);
        end
    end
end

function [ X_new ] = mutate_convlike( X )
%CROSSOVER_CONVLIKE 此处显示有关此函数的摘要
%   convolution-like operation

    [m, n] = size(X);

    X = padarray(X, [1 1], 'symmetric', 'both');
    tensor_X = twist(X);
    % block vectorize: bv_X=bvec(tensor_X);
 
   %% construct variation tensor of fixed size 
    l1 = m + 2;l2 = m + 2; l3 = n + 2; 
    A = zeros([l1, l2, l3]); % Creates a l1 x l2 x l3 tensor of zeros
    % Mean kernel
    for i = 2 : l1-1
        A(i, i-1:i+1, 1) = [1 1 1];
        A(i, i-1:i+1, 2) = [1 1 1];
    end
    A(:, :, l3) = A(:, :, 2);
    A = A / 9;
    % Guassian kernel
%     for i = 2 : l1-1
%         A(i, i-1:i+1, 1)=[2 4 2];
%         A(i, i-1:i+1, 2)=[1 2 1];
%     end
%     A(:, :, l3)=A(:, :, 2);
%     A = A / 16;
    X_new = tprod(A, tensor_X);
    X_new = squeeze(X_new);
    X_new = X_new(2 : end-1, 2 : end-1);
end