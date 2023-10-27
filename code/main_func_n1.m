
function [cost_best,pop_per_gen] = main_func_n1(param,F,C)

    popsize=param(:,1); % population size P
    Gen=param(:,2); % generation size
    R=param(:,3); % the percentage of replication of well-performed chromosomes  5%
    pc=param(:,4); % crossover rate
    pml=param(:,5); % mutation rate local
    pmg=param(:,6); % mutation rate global

    M=24; % number of machines
    m_grids=5; n_grids=6; % grid size

    pop=init_gen(popsize,m_grids,n_grids,M); % 初始化第一代
    fit= calc_fit(pop,F,C,M);

    pop_new=[]; %存储新的一代
    pop_per_gen=zeros(Gen+1,m_grids,n_grids,popsize);
    pop_per_gen(1,:,:,:)=pop;
    
    for t=1:Gen
        % 复制少量最优个体替换最差个体
        [pop, fit] = re_gen(pop,fit,R);
        %%
        % selection operator: select #popsize chromosomes from population
        pop=select_rssr(pop,fit);
        pop_new=crossover_n1(pop,pc); %crossover 
        pop_new=mutate_local(pop_new,pml); %mutate
        pop_new=mutate_global(pop_new,pmg); 

        %% update to the next generation
        fit_new=calc_fit(pop_new,F,C,M);
        pop=pop_new;fit=fit_new;
        pop_per_gen(t+1,:,:,:)=pop;
    end
    % final
    [~,best_idx]=max(fit_new);
    pop_best=pop_new(:,:,best_idx);
    cost_best = calc_cost(pop_best,F,C,M);
end

function [ pop ] = init_gen( popsize,m_grid,n_grid,M )
%INIT_GEN 此处显示有关此函数的摘要
%   此处显示详细说明
    pop=zeros(m_grid,n_grid,popsize);
    % num of objects to be processed
    nobjs=M; % if nobjs=m_grids*n_grids,判定是否存在空节点 -1表示一个restricted area，0表示 a dummy machine

    for i=1:popsize
        chrom=pop(:,:,i);
        for k=1:nobjs %从1到nobjs
            % 随机生成位置 (x,y),[1-m_grids],[1-n_grids]随机正整数
            x=unidrnd(m_grid);
            y=unidrnd(n_grid);
            while chrom(x,y)~=0 %判断该位置是否为空,禁止区或者被分配
                x=unidrnd(m_grid);
                y=unidrnd(n_grid);
            end
            chrom(x,y)=k;
        end
        pop(:,:,i)=chrom;
    end
end

function [ pop_new ] = crossover_n1( pop,pc )
    %crossover  交叉区域为行和列的混合区域，一次张量积
    [m,n,popsize]=size(pop); 
    pop_new=zeros(m,n,popsize);
     
    for k=1:2:(popsize-1)
        cid=randperm(popsize,2);
        % 两个染色体的三阶张量
        X(:,:,1)=pop(:,:,cid(1,1));
        X(:,:,2)=pop(:,:,cid(1,2));
        if rand<pc
            %变换张量,随机选择交换的行，列，会产生两个同时交换的区域
            %strating point of selected rows & columns
            sp_r=unidrnd(m); 
            sp_c=unidrnd(n);
            %number of selected rows & columns [1,m]&[1,n]
            L_r=unidrnd(m); L_c=unidrnd(n);
            diag_r=ones(m,1);diag_c=ones(n,1);
              
            for i=1:L_r
                if mod(sp_r+i-1,m)==0
                    diag_r(5,:)=0;
                else
                    diag_r(mod(sp_r+i-1,m),:)=0;
                end
            end
            
            for i=1:L_c
                if mod(sp_c+i-1,n)==0
                    diag_c(5,:)=0;
                else
                    diag_c(mod(sp_c+i-1,n),:)=0;
                end
            end
            
            A(:,:,1)=diag(diag_r);A(:,:,2)=diag(1-diag_r);
            B(:,:,1)=diag(diag_c);B(:,:,2)=diag(1-diag_c);
            
            if rand<0.5
                X_new=tprod(A,X);
            else
                X_new=tprod(X,B);
            end
           
           %% repair mechanism
            X_new(:,:,1)=repair_c(X_new(:,:,1),X(:,:,1));
            X_new(:,:,2)=repair_c(X_new(:,:,2),X(:,:,2));

            pop_new(:,:,k)=X_new(:,:,1);pop_new(:,:,k+1)=X_new(:,:,2);
        else
            pop_new(:,:,k)=X(:,:,1);pop_new(:,:,k+1)=X(:,:,2);
        end
    end
    
end

function [ pop_new , fit_new] = re_gen(pop,fitness,R)
%REPRODUCE 此处显示有关此函数的摘要
% remove the worst chromosomes floor(popsize*R)
% duplicating and inserting the best chromosomes into the current
% population

    [~,~,popsize]=size(pop);
    num=floor(popsize*R);
    pop_new=pop;
    fit_new=fitness;
    
    [~,ind_worst]=sort(fitness,'ascend'); % worst chromosomes
    [~,ind_best]=sort(fitness,'descend'); % best chromosomes
    
    pop_new(:,:,ind_worst(1:num,:))=pop(:,:,ind_best(1:num,:)); % 复制最优的个体替换最差的个体
    fit_new(ind_worst(1:num,:),:)=fitness(ind_best(1:num,:),:);

end

function [ pop_new ] = select_rssr(pop,fitness)
%SELECT 此处显示有关此函数的摘要
%   remainder stochastic sampling with replacement 无回放余数随机选择
% 计算群体中每个个体在下一代群体中的生存期望数目，整数部分为每个个体在下一代群体中的生存数目
    
    [m,n,popsize]=size(pop);
    
    %计算群体中每个个体在下一代群体中的生存期望数目    
%     num_expectation = zeros(popsize,1);
    num_expectation=popsize.*(fitness./sum(fitness));

    %取生存期望数目整数部分作为对应个体在下一代群体中的生存数目
    %sum(num_expectation_int)可以确定出下一代群体的个数
    num_expectation_int = floor(num_expectation);
    
    %用于存放新的个体
    pop_new=zeros(m,n,popsize); 
    
    %以下方法为按照生存期望选择的个体
    s=0;%计数，用于标记新的群体中个数
    for i=1:popsize
        if num_expectation_int(i)~=0
            a=[];%用于标识群体中的个体
            for j=1:num_expectation_int(i)
                if j==1
                    a(:,:,1)=pop(:,:,i);
                else
                    a(:,:,end+1)=pop(:,:,i);
                end
            end
            pop_new(:,:,s+1:s+num_expectation_int(i))=a;
            s=s+num_expectation_int(i);
        end
    end
    
    %计算新的适应度
    %确定下一代中还未确定的的 n - sum(num_expectation_int)个个体
    fitness_new=num_expectation-num_expectation_int; %小数部分用来计算轮盘赌的选择概率
%         fitness_new=fitness-sum(fitness).*num_expectation_int./popsize;
%     for i=1:popsize
%         fitness_new(i)=fitness(i) - num_expectation_int(i)*sum(fitness)/popsize; %num_expectation=popsize.*(fitness./sum(fitness))
%     end
    
    %以下为按照轮盘赌选择剩余的个体
    %以fitness_new为各个个体的新的适应度，再用基本的比例选择方法来随机确定下一代群体中剩余未确定的个体
    %Tips:当V为全零向量时，该算法无效，将随机选择个体;否则算法将从重要性指标不为0的个体中选择。
    %index:选的num_remain个个体的位置索引 
    
    %剩下的需要轮盘赌的个数 num_remain = popsize - s;
    num_remain=popsize-s;
    fitness_table = cumsum(fitness_new./sum(fitness_new));%累加 fitness_table = zeros(popsize,1);
    pop_remain=[];
    rs=sort(rand(num_remain,1));
    fiti=1;
    newi=1;
    while newi<=num_remain
        if rs(newi) < fitness_table(fiti)
            pop_remain(:,:,newi)=pop(:,:,fiti);
            newi=newi+1;
        else
            fiti=fiti+1;
        end
    end
    pop_new(:,:,s+1:end)=pop_remain;
end

function [ pop_new ] = mutate_local(pop,pm)
    [m,n,popsize]=size(pop); 
    pop_new=zeros(m,n,popsize);
        
    for i=1:popsize 
%         ind=unidrnd(popsize); %选择的父个体
        ind=i;
        if rand<pm %判断是否变异
            off=abs(round(mutate_convlike(pop(:,:,ind))));
           %% repair mechanism
            off=repair_m(off,pop(:,:,ind));
            pop_new(:,:,i)=off;
        else
            pop_new(:,:,i)=pop(:,:,ind); %insert the new chroms into a new population until the new population is full
        end
    end
end

function [ pop_new ] = mutate_global(pop,pm)
%MUTATE 单点变异
    ub=0;lb=24;
    [m,n,popsize]=size(pop); 
    mu_ind=rand(m,n,popsize)<pm;
    %mutant tensor
    gaussian_rand=normrnd(0,0.1,[m,n,popsize]).*(ub-lb);
    MT=mu_ind.*gaussian_rand;
    
    pop_new=abs(round(pop+MT));
    ind=[];zerom=zeros(m,n);
    for i=1:popsize
        if sum(sum(MT(:,:,i)))~=0
            ind=[ind;i];
        end
    end
    %repair
    for i=1:size(ind,1)
        pop_new(:,:,ind(i,:))=repair_m(pop_new(:,:,ind(i,:)),pop(:,:,ind(i,:)));
    end
end

function [ X_new ] = mutate_convlike(X)
%CROSSOVER_CONVLIKE 此处显示有关此函数的摘要
%   convolution-like operation
% 利用张量乘积的表示形式，构造的张量A根据问题定义设计

    [m,n]=size(X);%X大小m*n

    % 为了匹配张量运算的维度，将每个个体进行padding，四周填充0
    X=padarray(X,[1 1],'symmetric','both'); %大小变为m+2*n+2

    % 个体转换成tensor,即二维矩阵通过twist旋转转化成tensor 
    tensor_X=twist(X); %大小m*1*n

    % block vectorize: bv_X=bvec(tensor_X);
 
   %% 构造transform tensor，大小l1,l2,l3
    l1=m+2;l2=m+2;l3=n+2; 
    A=zeros(l1,l2,l3); % Creates a l1 x l2 x l3 tensor of zeros
    %均值核
    for i=2:l1-1
        A(i,i-1:i+1,1)=[1 1 1];
        A(i,i-1:i+1,2)=[1 1 1];
    end
    A(:,:,l3)=A(:,:,2);
    A=A/9;
    %高斯核
%     for i=2:l1-1
%         A(i,i-1:i+1,1)=[2 4 2];
%         A(i,i-1:i+1,2)=[1 2 1];
%     end
%     A(:,:,l3)=A(:,:,2);
%     A=A/16;
    % block circutant matrix: bc_A=bcirc(A);
    %t-product
    X_new=tprod(A,tensor_X);

    % 将新得到的tensor旋转得到矩阵，然后从中提取unpadding个体
    X_new=squeeze(X_new);
    X_new=X_new(2:end-1,2:end-1);

end

% 高斯（0,1) or (0,0.1)
