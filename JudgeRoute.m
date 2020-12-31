%
%      @作者：随心390
%      @微信公众号：优化算法交流地
%
%% 判断当前路径是否满足时间窗约束和载重量约束，0表示违反约束，1表示满足全部约束
%输入：chrom               个体
%输入：cap                 最大载重量
%输入：demands             需求量
%输入：a                   顾客时间窗开始时间[a[i],b[i]]
%输入：b                   顾客时间窗结束时间[a[i],b[i]]
%输入：L                   配送中心时间窗结束时间
%输入：s                   客户点的服务时间
%输入：dist                距离矩阵，满足三角关系，暂用距离表示花费c[i][j]=dist[i][j]
%输出：flag                0表示违反约束，1表示满足全部约束
function [flag,Q_U,wait,w_pp]=JudgeRoute(route,cap,demands,a,b,L,s,dist)
flag=1;                         %假设满足约束
lr=length(route);               %该条路径上顾客数目
%% 计算每辆车的装载量
Ld=leave_load(route,demands);
%% 计算载重率
% p_b 初始油耗
p_b=0.18;
P_e=0.41;
% p_e 满载油耗
Q_U=p_b+Ld*(P_e-p_b)/cap; %耗油能力 
wait=0;
w_pp=1;
%如果满足载重量约束，需用进行时间窗判断
if Ld<=cap
    %% 计算该路径上在各个点开始服务的时间，还计算返回配送中心时间
    [arr,bs,wait,back]=begin_s(route,a,s,dist);
    
    if wait(1)~=0
        w_pp=wait(1)/(b(route(1))-bs(1));
    else
        w_pp=1;
    end
    %如果满足配送中心右时间窗约束，需用进行判断各个顾客的时间窗是否满足时间窗约束
    if back<=L
        for i=1:lr
            
            %一旦发现某个顾客的时间窗是否满足时间窗约束，则直接判为违反约束，将flag设为0
            
            if bs(i)>b(route(i))
                flag=0;
            end
           if wait(i)~=0
                w_pp(i)=(b(route(i))-bs(i))/wait(i);
            else
                w_pp(i)=1;
            end
        end
    else
        %如果不满足配送中心右时间窗约束，直接判为违反约束，将flag设为0
        flag=0;
    end
else
    %如果不满足载重量约束，不用进行时间窗判断，直接判为违反约束，将flag设为0
    flag=0;
end
end