
function [ tensor_X ] = twist( X )
%TWIST �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    % �������X��Ӧ������
    % ���룺����X����Сm*n
    % ���������tensor_X, ��Сm*1*n
    [m,n]=size(X);
    tensor_X = reshape(X,m,1,n);

end

