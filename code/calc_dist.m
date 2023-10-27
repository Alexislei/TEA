function [ dist ] = calc_dist(chrom,M)
%CALC_DIST 此处显示有关此函数的摘要
%   此处显示详细说明
    %从1到M表示每台机器
    dist=zeros(M,M); %M:number of machines
    % rectangular distance 
    for i=1:M
        [x1,y1]=find(chrom==i);
        if isempty(x1)
            dist(i,:)=0;
        else
            for j=1:M
                [x2,y2]=find(chrom==j);
                if isempty(x2)
                    dist(i,j)=0;
                else
                    dist(i,j)=abs(x1-x2)+abs(y1-y2);
                end
            end
        end
    end  
end

