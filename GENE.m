function [final_solution,gbest]=GENE(pop,maxite,num_cargo,pm,box)
n=pop;
ger=maxite;
dim=num_cargo;
pm;
% n-- 种群规模% ger-- 最大迭代次数%pm-- 变异概率
% v-- 初始种群（规模为n）% fit-- 适应度向量


global num_box;


% 生成初始种群
for i=1:n
    v(i,:)= fix((num_box)*rand(1,dim))+1; 
    Scheme=transform(v(i,:));                     %解转化成“货箱：货物”对应的形式
    [feas_solution,Scheme]= placement(Scheme,box);             %装箱处理
    [PG,PV,pbest ]= evaluate(feas_solution) ;      %计算适应度
    fit(i)=pbest;
end

[gbest, index] = max(fit);
    for i=1:n
        v2(i,:) = v(index,:);
    end
    
[N,L]=size(v);           %得到初始规模行，列

it=1; % 迭代计数器
% 开始进化
while it<=ger 
    
    %选择最佳保留复制    
    for i=1:n
        if i~=index
            if rand > pm    % 交叉
                begin = 1 + fix(rand*dim);
                stop = begin + fix((dim-begin)*rand);
                v2(i,begin:stop) = v(i,begin:stop);
                
            else            % 变异
                begin = 1 + fix(rand*dim);
                stop = begin + fix((dim-begin)*rand);
                v2(i,begin:stop)=fix(rand(1,(stop-begin+1))*num_box)+1;
            end
        else
            continue;
        end
         final_solution=v2(1,:);
        Scheme=transform(v2(i,:));                   
        [feas_solution,Scheme]= placement(Scheme,box);             
        [PG,PV,fit(i) ]= evaluate(feas_solution) ;      
        if fit(i) > gbest
            gbest = fit(i);
            index = i;
            final_solution = v2(i,:);
        else
           
        end
            
    end
    
    v=v2;   %复制
    for i=1:n
        v2(i,:) = v(index,:);
    end
    
    it=it+1;            
end      

