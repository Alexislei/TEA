function [ fitness ] = calc_fit( pop,F,C,M)
%
    [~,~,popsize]=size(pop);
    fitness=zeros(popsize,1); % 适应值保存为popsize*1的向量
    for i=1:popsize
%         pop(:,:,i)
        D=calc_dist(pop(:,:,i),M);% rectangular distance
        totalCost=sum(sum(F.*C.*D));
        fitness(i,:)=totalCost;
    end
    fitness=1./fitness; % 错误判断：是否出现fitness=0的情况
end

