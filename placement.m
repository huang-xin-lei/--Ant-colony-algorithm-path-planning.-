function [ newsolution,scheme ] = placement( scheme,box )

global cargo; global num_box;global solution
newsolution=solution;
num=size(scheme,2);
for i=1:num_box
    x=0;y=0;z=0;   % 新货装载位置起点
    xb=0;yb=0;zb=0;   % 当前边界
    j=1;
    while (j <= num )
      if scheme(i,j)~=0
          if cargo(scheme(i,j),4)>box(i,4) || cargo(scheme(i,j),3)>box(i,3) ||cargo(scheme(i,j),2)>box(i,2)
               newsolution( scheme(i,j))=0;
               scheme(i,j)=0;
               j=j+1;
          else
              if (cargo(scheme(i,j),4)+z ) < box(i,4)   % Z越界  则判定已装不下
                if (cargo(scheme(i,j),3)+y ) < box(i,3)  %Z未越界 Y越界 则考虑上一层
                    if (cargo(scheme(i,j),2)+x ) < box(i,2) %X Y Z 均未越界，当前货物装载成功，否则考虑下一排
                        x = x + cargo(scheme(i,j),2);                     
                        nyb = y + cargo(scheme(i,j),3);
                        nzb = z + cargo(scheme(i,j),4);
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
               newsolution( scheme(i,j))=0;
               scheme(i,j)=0;
               j=j+1;
              end
          end
      else
         j=j+1; 
      end          
  end
end


