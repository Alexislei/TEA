
function [cost_best,pop_per_gen] = main_func_n1(param,F,C)

    popsize=param(:,1); % population size P
    Gen=param(:,2); % generation size
    R=param(:,3); % the percentage of replication of well-performed chromosomes  5%
    pc=param(:,4); % crossover rate
    pml=param(:,5); % mutation rate local
    pmg=param(:,6); % mutation rate global

    M=24; % number of machines
    m_grids=5; n_grids=6; % grid size

    pop=init_gen(popsize,m_grids,n_grids,M); % ��ʼ����һ��
    fit= calc_fit(pop,F,C,M);

    pop_new=[]; %�洢�µ�һ��
    pop_per_gen=zeros(Gen+1,m_grids,n_grids,popsize);
    pop_per_gen(1,:,:,:)=pop;
    
    for t=1:Gen
        % �����������Ÿ����滻������
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
%INIT_GEN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    pop=zeros(m_grid,n_grid,popsize);
    % num of objects to be processed
    nobjs=M; % if nobjs=m_grids*n_grids,�ж��Ƿ���ڿսڵ� -1��ʾһ��restricted area��0��ʾ a dummy machine

    for i=1:popsize
        chrom=pop(:,:,i);
        for k=1:nobjs %��1��nobjs
            % �������λ�� (x,y),[1-m_grids],[1-n_grids]���������
            x=unidrnd(m_grid);
            y=unidrnd(n_grid);
            while chrom(x,y)~=0 %�жϸ�λ���Ƿ�Ϊ��,��ֹ�����߱�����
                x=unidrnd(m_grid);
                y=unidrnd(n_grid);
            end
            chrom(x,y)=k;
        end
        pop(:,:,i)=chrom;
    end
end

function [ pop_new ] = crossover_n1( pop,pc )
    %crossover  ��������Ϊ�к��еĻ������һ��������
    [m,n,popsize]=size(pop); 
    pop_new=zeros(m,n,popsize);
     
    for k=1:2:(popsize-1)
        cid=randperm(popsize,2);
        % ����Ⱦɫ�����������
        X(:,:,1)=pop(:,:,cid(1,1));
        X(:,:,2)=pop(:,:,cid(1,2));
        if rand<pc
            %�任����,���ѡ�񽻻����У��У����������ͬʱ����������
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
%REPRODUCE �˴���ʾ�йش˺�����ժҪ
% remove the worst chromosomes floor(popsize*R)
% duplicating and inserting the best chromosomes into the current
% population

    [~,~,popsize]=size(pop);
    num=floor(popsize*R);
    pop_new=pop;
    fit_new=fitness;
    
    [~,ind_worst]=sort(fitness,'ascend'); % worst chromosomes
    [~,ind_best]=sort(fitness,'descend'); % best chromosomes
    
    pop_new(:,:,ind_worst(1:num,:))=pop(:,:,ind_best(1:num,:)); % �������ŵĸ����滻���ĸ���
    fit_new(ind_worst(1:num,:),:)=fitness(ind_best(1:num,:),:);

end

function [ pop_new ] = select_rssr(pop,fitness)
%SELECT �˴���ʾ�йش˺�����ժҪ
%   remainder stochastic sampling with replacement �޻ط��������ѡ��
% ����Ⱥ����ÿ����������һ��Ⱥ���е�����������Ŀ����������Ϊÿ����������һ��Ⱥ���е�������Ŀ
    
    [m,n,popsize]=size(pop);
    
    %����Ⱥ����ÿ����������һ��Ⱥ���е�����������Ŀ    
%     num_expectation = zeros(popsize,1);
    num_expectation=popsize.*(fitness./sum(fitness));

    %ȡ����������Ŀ����������Ϊ��Ӧ��������һ��Ⱥ���е�������Ŀ
    %sum(num_expectation_int)����ȷ������һ��Ⱥ��ĸ���
    num_expectation_int = floor(num_expectation);
    
    %���ڴ���µĸ���
    pop_new=zeros(m,n,popsize); 
    
    %���·���Ϊ������������ѡ��ĸ���
    s=0;%���������ڱ���µ�Ⱥ���и���
    for i=1:popsize
        if num_expectation_int(i)~=0
            a=[];%���ڱ�ʶȺ���еĸ���
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
    
    %�����µ���Ӧ��
    %ȷ����һ���л�δȷ���ĵ� n - sum(num_expectation_int)������
    fitness_new=num_expectation-num_expectation_int; %С�����������������̶ĵ�ѡ�����
%         fitness_new=fitness-sum(fitness).*num_expectation_int./popsize;
%     for i=1:popsize
%         fitness_new(i)=fitness(i) - num_expectation_int(i)*sum(fitness)/popsize; %num_expectation=popsize.*(fitness./sum(fitness))
%     end
    
    %����Ϊ�������̶�ѡ��ʣ��ĸ���
    %��fitness_newΪ����������µ���Ӧ�ȣ����û����ı���ѡ�񷽷������ȷ����һ��Ⱥ����ʣ��δȷ���ĸ���
    %Tips:��VΪȫ������ʱ�����㷨��Ч�������ѡ�����;�����㷨������Ҫ��ָ�겻Ϊ0�ĸ�����ѡ��
    %index:ѡ��num_remain�������λ������ 
    
    %ʣ�µ���Ҫ���̶ĵĸ��� num_remain = popsize - s;
    num_remain=popsize-s;
    fitness_table = cumsum(fitness_new./sum(fitness_new));%�ۼ� fitness_table = zeros(popsize,1);
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
%         ind=unidrnd(popsize); %ѡ��ĸ�����
        ind=i;
        if rand<pm %�ж��Ƿ����
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
%MUTATE �������
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
%CROSSOVER_CONVLIKE �˴���ʾ�йش˺�����ժҪ
%   convolution-like operation
% ���������˻��ı�ʾ��ʽ�����������A�������ⶨ�����

    [m,n]=size(X);%X��Сm*n

    % Ϊ��ƥ�����������ά�ȣ���ÿ���������padding���������0
    X=padarray(X,[1 1],'symmetric','both'); %��С��Ϊm+2*n+2

    % ����ת����tensor,����ά����ͨ��twist��תת����tensor 
    tensor_X=twist(X); %��Сm*1*n

    % block vectorize: bv_X=bvec(tensor_X);
 
   %% ����transform tensor����Сl1,l2,l3
    l1=m+2;l2=m+2;l3=n+2; 
    A=zeros(l1,l2,l3); % Creates a l1 x l2 x l3 tensor of zeros
    %��ֵ��
    for i=2:l1-1
        A(i,i-1:i+1,1)=[1 1 1];
        A(i,i-1:i+1,2)=[1 1 1];
    end
    A(:,:,l3)=A(:,:,2);
    A=A/9;
    %��˹��
%     for i=2:l1-1
%         A(i,i-1:i+1,1)=[2 4 2];
%         A(i,i-1:i+1,2)=[1 2 1];
%     end
%     A(:,:,l3)=A(:,:,2);
%     A=A/16;
    % block circutant matrix: bc_A=bcirc(A);
    %t-product
    X_new=tprod(A,tensor_X);

    % ���µõ���tensor��ת�õ�����Ȼ�������ȡunpadding����
    X_new=squeeze(X_new);
    X_new=X_new(2:end-1,2:end-1);

end

% ��˹��0,1) or (0,0.1)
