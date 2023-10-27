
function [ tensor_X ] = twist( X )
%TWIST 此处显示有关此函数的摘要
%   此处显示详细说明
    % 构造矩阵X对应的张量
    % 输入：矩阵X，大小m*n
    % 输出：张量tensor_X, 大小m*1*n
    [m,n]=size(X);
    tensor_X = reshape(X,m,1,n);

end

