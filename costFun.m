%
%      @作者：随心390
%      @微信公众号：优化算法交流地
%
%% 计算一个配送方案的总成本=cd*车辆使用数目+ct*车辆行驶总距离
%输入VC：          配送方案
%输入dist：        距离矩阵
%输出cost：        总成本
%输出NV：          车辆使用数目
function [cost,NV,TD,ff,c_11,w_ppm,w_pp]=costFun(VC,dist,Q_U,w_pp)
NV=size(VC,1);                      %车辆使用数目
TD=travel_distance(VC,dist);        %行驶总距离
ff=mean(Q_U)*TD;%耗油量
u=30;%碳税
w=2.669;%碳排放系数
cost=150*NV+NV*ff*5.41+NV*u*w*ff;%固定成本+油耗成本+碳税成本；150 单车固定成本；5.41 油价
c_11=NV*w*ff; %碳排放量
w_ppm=mean(w_pp);%平均满意度
if w_ppm>1
    w_ppm=1;
end
end