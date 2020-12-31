%% 源代码参考 
% 网上来源 具体链接后续补上
% 修改 尚帝
% 内容： 将时间窗口约束的路径规划结合背包装箱问题，将装箱问题可视化。
% 特点： 在优化过程中考虑背包的装置也摆放位置，多角度考虑问题，细化过程和可视化
% 不足： 没有真正对解决问题提供参考意义，缺少实际价值。没有将两者更加有机结合
%待改进的点： 1.将装箱部分与路径选择的权重相关联，考虑摆放位置与用户体验，运送货物的需求相匹配。可视化的过程与具体要求匹配
% 2. 降低参数设置，简化算法
% 邮箱 162480875@qq.com 有疑问或建议发送

%% 2020年12月19日16:34:43
% 版本 0.01
% 加入了满意度 碳排放量
%% 2020年12月27日15:51:32
% 版本 0.02
% 简单调整了作图增加了标签
%% 2020年12月30日23:11:18
% 版本0.1
% 增加了装箱过程展示与装箱优化和车辆类型选择
%% 2020年12月31日10:24
% 版本0.11
% 将集装箱约束与车辆约束嵌套
clear
clc
close all
tic
%% 用importdata这个函数来读取文件
c101=importdata('c103.txt');
box=load('box'); %获箱体信息  box: 重 长 宽 高 体积
box_sty=1; %装载箱类型种类
orginal_cargo=load('cargo'); %获取不同地点的货物数目与货物属性  序号 货物数目 重 长 宽 高
cap=box(1);                                                        %车辆最大装载量
v_cap=2;   %车辆速度（单位时间）
%% 提取数据信息
E=c101(1,5);                                                    %配送中心时间窗开始时间
L=c101(1,6);                                                    %配送中心时间窗结束时间
vertexs=c101(:,2:3);                                            %所有点的坐标x和y
customer=vertexs(2:end,:);                                      %顾客坐标
cusnum=size(customer,1);                                        %顾客数
v_num=3;                                                       %车辆最多使用数目
demands=c101(2:end,4);                                          %需求量
a=c101(2:end,5);                                                %顾客时间窗开始时间[a[i],b[i]]
b=c101(2:end,6);                                                %顾客时间窗结束时间[a[i],b[i]]
width=b-a;                                                      %顾客的时间窗宽度
s=c101(2:end,7);                                                %客户点的服务时间
choose_way=input('计算的距离类型————1 坐标位置 2 经纬度===');
if choose_way==1
    h=pdist(vertexs);
    dist=squareform(h);                                             %距离矩阵 坐标相减
   % dist(1,3)=555;
else 
     mi=[vertexs];
        cc=zeros(size(mi,1));
        for j=1:size(mi,1)
            for i=j:size(mi,1)
                cc(j,i)=abs(6371004*acos((sin(deg2rad(mi(j,2)))*sin(deg2rad(mi(i,2)))+cos(deg2rad(mi(j,2)))*cos(deg2rad(mi(i,2)))*cos(deg2rad(mi(i,1)-mi(j,1))))));
                if i==j || cc(j,i)==0
                     cc(j,i)=eps;
                end
                cc(i,j)=cc(j,i);
            end
        end
        un_1_1=cc;
        dist=un_1_1;        %距离矩阵
end
%% 初始化参数
w_PPm_b=1;%初始满意度为1
m=50;                                                           %蚂蚁数量
alpha=1;                                                        %信息素重要程度因子
beta=3;                                                         %启发函数重要程度因子
gama=2;                                                         %等待时间重要程度因子
delta=3;                                                        %时间窗跨度重要程度因子
r0=0.5;                                                         %r0为用来控制转移规则的参数
rho=0.85;                                                       %信息素挥发因子
Q=100;                                                          %更新信息素浓度的常数
Eta=1./dist;                                                    %启发函数
Tau=ones(cusnum+1,cusnum+1);                                    %信息素矩阵
Table=zeros(m,cusnum);                                          %路径记录表
iter=1;                                                         %迭代次数初值
iter_max=10;                                                   %最大迭代次数
Route_best=zeros(iter_max,cusnum);                              %各代最佳路径
Cost_best=zeros(iter_max,1);                                    %各代最佳路径的成本
ture_c=0;  %是否满足装箱 
%% 迭代寻找最佳路径
while iter<=iter_max || ture_c~=1
    %% 先构建出所有蚂蚁的路径
    %逐个蚂蚁选择
    for i=1:m
        %逐个顾客选择
        for j=1:cusnum
            r=rand;                                             %r为在[0,1]上的随机变量
            np=next_point(i,Table,Tau,Eta,alpha,beta,gama,delta,r,r0,a,b,width,s,L,dist,cap,demands);
            Table(i,j)=np;
        end
    end
    %% 计算各个蚂蚁的成本=1000*车辆使用数目+车辆行驶总距离
    cost=zeros(m,1);
    NV=zeros(m,1);
    TD=zeros(m,1);
    for i=1:m
        [VC,NV,TD,Q_U,w_pp]=decode(Table(i,:),cap,demands,a,b,L,s,dist);
        [cost(i,1),NV(i,1),TD(i,1),ff(i,1),c_11(i,1),w_PPm(i,1)]=costFun(VC,dist,Q_U,w_pp);
    end
    %% 计算最小成本及平均成本
    if iter == 1
        [min_Cost,min_index]=min(cost);
        c_bb=c_11(min_index);
        f_bb=ff(min_index);
        Cost_best(iter)=min_Cost;
        w_PPm_m(iter)=mean(w_PPm);
        Route_best(iter,:)=Table(min_index,:);
    else
        [min_Cost,min_index]=min(cost);
        w_PPm_m(iter)=mean(w_PPm);
        Cost_best(iter)=min(Cost_best(iter - 1),min_Cost);
        c_bb=c_11(min_index);
        f_bb=ff(min_index);
        w_PPm_b=w_PPm(min_index);
        if Cost_best(iter)==min_Cost
            Route_best(iter,:)=Table(min_index,:);
            c_bb=c_11(min_index);
            f_bb=ff(min_index);
             w_PPm_b=w_PPm(min_index);
        else
            Route_best(iter,:)=Route_best((iter-1),:);
            c_bb=c_11(min_index);
            f_bb=ff(min_index);
             w_PPm_b=w_PPm(min_index);
        end
    end
    %% 更新信息素
    bestR=Route_best(iter,:);
    [bestVC,bestNV,bestTD]=decode(bestR,cap,demands,a,b,L,s,dist);
    Tau=updateTau(Tau,bestR,rho,Q,cap,demands,a,b,L,s,dist);
   %% 判断生成的路径是否符合要求
   % 测试功能
    for i=1:size(bestVC,1)
         aa=bestVC{i,:};
        [PG,PV,gbest,timecost,Scheme]=Container(aa,box,orginal_cargo,box_sty); %对装修进行约束
            if PG>1||PV>1 % 判断满载率和空间载重限制情况 如果不满足转化货运箱类型
                box_sty=box_sty+1;
                ture_c=0;
                 disp('没有符合要求的装箱类型,需要改变种类')
%                  iter=abs(iter-1);
                kk=1;
                if max(box_sty)>size(box,1) 
                     disp('装箱类型全部都不符合要求,重新设定：已暂停')
%                     pause
                      box_sty=[1:kk];
                      kk=kk+1;
                      ture_c=0;
                end
            else
                ture_c=1;
            end
    end
    %% 最短时间
    t_1=bestTD/v_cap; % 路径时间
    s_1=sum(s); % 总服务时间
    T_BEAST=t_1+s_1; 
    %% 打印当前最优解
    if v_num<=num2str(bestNV)
        disp(['第',num2str(iter),'代最优解:'])
        disp(['车辆使用数目：',num2str(bestNV),'，装载箱种类选择：',num2str(box_sty),'，车辆行驶总距离：',num2str(bestTD),'最低碳排放量',num2str(c_bb),'满意度为',num2str(w_PPm_b),'最短时间',num2str(T_BEAST)]);
    else
        disp(['第',num2str(iter),'代最优解:'])
        disp(['车辆循环次数：',num2str(bestNV),'，装载箱种类选择：',num2str(box_sty),'，车辆行驶总距离：',num2str(bestTD),'最低碳排放量',num2str(c_bb),'满意度为',num2str(w_PPm_b),'最短时间',num2str(T_BEAST)]);
    end
    fprintf('\n')
   
   %% 迭代次数加1，清空路径记录表
    iter=iter+1;
    Table=zeros(m,cusnum);
end
%% 结果显示
bestRoute=Route_best(end,:);
[bestVC,NV,TD]=decode(bestRoute,cap,demands,a,b,L,s,dist);
draw_Best(bestVC,vertexs);
 for i=1:size(bestVC,1)
%     currentDate=strcat('a',num2str(i))
%     aa.(currentDate)=bestVC{i,:};
        aa=bestVC{i,:};
%     save(strcat('bestVC',num2str(i),'.mat'), 'aa')
[PG,PV,gbest,timecost,Scheme]=Container(aa,box,orginal_cargo,box_sty); %对装修进行约束
result(Scheme,15);      %将装箱方案Scheme 按每行15个货物显示
figure
fprintf('重量利用率：\t%5.3f %%\n',PG*100);
fprintf('空间利用率：\t%5.3f %%\n',PV*100);
fprintf('综合利用率：\t%5.3f %%\n',gbest*100);
% fprintf('计算时间：\t\t%5.4f s\n',timecost);
% disp('图像生成中...')

        depict( Scheme, 1,'r' )   %    ( 方案，画出编号为i箱子，颜色） 颜色：r\g\b\c\m\y\k\w
 end

%% 绘图
figure
plot(1:iter_max,Cost_best(1:iter_max),'b')
xlabel('迭代次数')
ylabel('成本')
title('各代最小成本变化趋势图')
% %% 判断最优解是否满足时间窗约束和载重量约束，0表示违反约束，1表示满足全部约束
 flag=Judge(bestVC,cap,demands,a,b,L,s,dist);
%% 检查最优解中是否存在元素丢失的情况，丢失元素，如果没有则为空
DEL=Judge_Del(bestVC);
figure
plot(1:iter_max,w_PPm_m(1:iter_max),'b')
xlabel('迭代次数')
ylabel('满意度')
title('各代满意度平均变化趋势图')
figure
plot(1:m,w_PPm(1:m),'b')
xlabel('蚂蚁')
ylabel('满意度')
title('最后一代满意度平均变化趋势图')
toc