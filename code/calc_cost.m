
function [ cost ] = calc_cost(pop,F,C,M)
%CALC_COST 此处显示有关此函数的摘要

    [~,~,popsize]=size(pop);
    cost=zeros(popsize,1); % 适应值保存为popsize*1的向量
    for i=1:popsize
%         pop(:,:,i)
        D=calc_dist(pop(:,:,i),M); % rectangular distance
        totalCost=sum(sum(F.*C.*D));
        cost(i,:)=totalCost;
    end
end

