function [PG,PV,gbest,timecost,Scheme]=Container(aa,box,orginal_cargo,box_sty)
 %% 进行集装箱优化
global cargo; global lambda; global num_cargo;global num_box;global solution;

%% -------------------------------控制参数---------------------------

lambda = 0.5;       % 重量利用率权重

T0 = 100;           % 初始温度
T_End = 1;          % 终止温度
metropolis = 100;   % 退火算法中 metropolis链长度
cooling = 0.98;     % 降温系数

pop = 11;           %遗传算法染色体数
maxite = 10;       %遗传最大迭代次数
pm = 0.1;           %遗传变异概率
%% --------------------------------------------------------------------

%% ----------------------------初始化：读取货箱信息 ----------------------------

% aa=load(bestway); %获取路径上的解
box=box(box_sty,:); % 选择的种类
% cargo=[];
orginal_cargo=orginal_cargo(aa,:); 
count=1;
for i=1:size(orginal_cargo,1)           %重构货物格式  cargo: 重 长 宽 高 体积 ；其中 长>宽>高
    for j=1:orginal_cargo(i,2) %货物个数 =size(orginal_cargo,1)*orginal_cargo(i,2) 
        cargo(count,1:4) = orginal_cargo(i,3:6);
        cargo(count,5) = prod(cargo(count,2:4),2); 
        cargo(count,2:4) = sort(cargo(count,2:4),'descend');
        count=count+1;
    end
end         
for i=1:size(box,1)                          %重构箱子box: 重 长 宽 高 体积
    box(i,5)=prod(box(i,2:4),2);            
end

num_cargo=size(cargo,1);  % 货物数
num_box=size(box,1);      % 货箱数

solution= fix((num_box)*rand(1,num_cargo))+1;   %随机生成初始解
Scheme=transform(solution);                     %解转化成“货箱：货物”对应的形式
[feas_solution,Scheme]= placement(Scheme,box);             %装箱处理

[PG,PV,gbest ]= evaluate(feas_solution,box) ;      %计算适应度

%--------------------------------------------------------------------

%----------------------------退火------------------------
begin=cputime;   %开始计时

%遗传算法优化     GENE（染色体数/种群规模，最大迭代次数，染色体长度/维度，变异概率）
[final_solution,gbest]=GENE(pop,maxite,num_cargo,pm,box) ;  

%遗传执行完毕后  模拟退火进一步优化
T = T0;
while T > T_End
    for i=1:metropolis
        %-----------随机交换两件货物生成新解
        newsolution=final_solution;
        R1=fix(rand*num_cargo)+1;
        R2=fix(rand*num_cargo)+1;
        inter=newsolution(R1);
        newsolution(R1)=newsolution(R2);
        newsolution(R2)=inter;
        NewScheme=transform(newsolution);                   % 分配货箱
        [feas_solution,NewScheme]= placement(NewScheme,box);              % 装箱处理
        [NPG,NPV,pbest ]= evaluate(feas_solution,box);            % 评估新方案
        if pbest>gbest
            gbest = pbest;
            final_solution = newsolution;
            PG = NPG;
            PV = NPV;
            Scheme = NewScheme;
        else
            if  rand < exp( (pbest-gbest)*100*T0/T)
                gbest=pbest;
                final_solution=newsolution;
                PG = NPG;
                PV = NPV;
                Scheme = NewScheme;
            end
        end   
    end
    T = T * cooling;
end

timecost = cputime-begin;   %计时结束



