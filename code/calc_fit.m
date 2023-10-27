function [ fitness ] = calc_fit( pop,F,C,M)
%
    [~,~,popsize]=size(pop);
    fitness=zeros(popsize,1); % ��Ӧֵ����Ϊpopsize*1������
    for i=1:popsize
%         pop(:,:,i)
        D=calc_dist(pop(:,:,i),M);% rectangular distance
        totalCost=sum(sum(F.*C.*D));
        fitness(i,:)=totalCost;
    end
    fitness=1./fitness; % �����жϣ��Ƿ����fitness=0�����
end

