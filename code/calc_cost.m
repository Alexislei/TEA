
function [ cost ] = calc_cost(pop,F,C,M)
%CALC_COST �˴���ʾ�йش˺�����ժҪ

    [~,~,popsize]=size(pop);
    cost=zeros(popsize,1); % ��Ӧֵ����Ϊpopsize*1������
    for i=1:popsize
%         pop(:,:,i)
        D=calc_dist(pop(:,:,i),M); % rectangular distance
        totalCost=sum(sum(F.*C.*D));
        cost(i,:)=totalCost;
    end
end

