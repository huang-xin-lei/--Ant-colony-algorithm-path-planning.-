function [  ] = depict( scheme, index ,mycolor)

global cargo;global num_box;global box;

num=size(scheme,2);
vert=zeros(8,3);
%------------定义画出函数画出每个面----------------
     fac = [1 2 3 4; ...
        2 6 7 3; ...
        4 3 7 8; ...
        1 5 8 4; ...
        1 2 6 5; ...
        5 6 7 8];
    vert(2,:)=[0 box(index,3) 0];
    vert(3,:)=[box(index,2) box(index,3) 0];
    vert(4,:)=[box(index,2) 0 0];
    vert(5,:)=[0 0 box(index,4)];
    vert(6,:)=[0 box(index,3) box(index,4)];
    vert(7,:)=[box(index,2) box(index,3) box(index,4)];
    vert(8,:)=[box(index,2) 0 box(index,4)];
    patch('Faces',fac,'Vertices',vert,'FaceColor',mycolor);  % patch function
    alpha('color');
    view(30,30);
    hold on;
    
    x=0;y=0;z=0;   % 新货装载位置起点
    xb=0;yb=0;zb=0;   % 当前边界
    j=1;
    while (j <= num )
      if scheme(index,j)~=0
          if (cargo(scheme(index,j),4)+z ) < box(index,4)   % Z越界  则判定已装不下
             if (cargo(scheme(index,j),3)+y ) < box(index,3)  %Z未越界 Y越界 则考虑上一层
                 if (cargo(scheme(index,j),2)+x ) < box(index,2) %X Y Z 均未越界，当前货物装载成功，否则考虑下一排
                     vert(1,:)=[x y z];
                     vert(2,:)=[x y+cargo(scheme(index,j),3) z];
                     vert(3,:)=[x+cargo(scheme(index,j),2) y+cargo(scheme(index,j),3) z];
                     vert(4,:)=[x+cargo(scheme(index,j),2) y z];
                     vert(5,:)=[x y z+cargo(scheme(index,j),4)];
                     vert(6,:)=[x y+cargo(scheme(index,j),3) z+cargo(scheme(index,j),4)];
                     vert(7,:)=[x+cargo(scheme(index,j),2) y+cargo(scheme(index,j),3) z+cargo(scheme(index,j),4)];
                     vert(8,:)=[x+cargo(scheme(index,j),2) y z+cargo(scheme(index,j),4)];
                     patch('Faces',fac,'Vertices',vert,'FaceColor',mycolor);  
                     material shiny;
                     alpha('color');
                     alphamap('rampdown');
                     x = x + cargo(scheme(index,j),2);  
                     nyb = y + cargo(scheme(index,j),3);
                     nzb = z + cargo(scheme(index,j),4);
                     j=j+1;   
                     if nyb > yb
                         yb=nyb;
                     end
                     if nzb > zb
                         zb=nzb;
                     end
                else
                   x=0;y=yb;    
                end
             else
                x=0;y=0;z=zb;   % 新货装载位置起点 (向上摞一层)
                xb=0;yb=0;
             end
          else
               newsolution( scheme(index,j))=0;
               scheme(index,j)=0;
               j=j+1;
          end

     else
        j=j+1; 
     end          
   end 
    
end

