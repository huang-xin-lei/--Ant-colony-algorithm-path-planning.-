
%% 配送方案与路径之间进行转换
function Route=change(VC)
NV=size(VC,1);                         %车辆使用数目
Route=[];
for i=1:NV
    Route=[Route,VC{i}];
end
end

