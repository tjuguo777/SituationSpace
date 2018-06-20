function [C,a_max] = XYZ_frame_change(u,ngz)
    thea=asin(u(3)/norm(u));
    if u(1)>=0
        dem=0;    
    else
        if u(2)>=0      %�ڶ�����Ҫ��pi
            dem=1;
        else
            dem=-1;     %��������Ҫ��pi
        end  
    end
    fai=atan(u(2)/u(1))+pi*dem;
    C=[sin(fai) -cos(fai)*cos(thea) -cos(fai)*sin(thea);
       cos(fai)  sin(fai)*cos(thea)  sin(fai)*sin(thea); 
       0        -sin(thea)                    cos(thea)];
    a_max=A_max_get(ngz,thea);
end

