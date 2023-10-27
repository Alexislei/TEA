
function [ fit ] = evaluate(pop,F,C,M)
% Evaluation function: evaluate the fitness of each individual
    [~,~,popsize]=size(pop);
    fit=zeros(popsize,1); % fitness vector
    for i=1:popsize 
        cost=calc_fit(pop(:,:,i),M); % define the fitness function based on the target problem
        fit(i,:)=cost;
    end
end

