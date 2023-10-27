
function [ tensor_X ] = twist( X )
% twist function: rotate a 2-dim individual to a 3-order tensor while preserving numerical property
% input: a two-dim individual，大小m*n
% output: a third-order individual tensor_X, size m*1*n
    [m,n]=size(X);
    tensor_X = reshape(X,m,1,n);
end

