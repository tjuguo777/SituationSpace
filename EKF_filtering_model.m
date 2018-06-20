function [X0,P0,v,Pzminus]=EKF_filtering_model(F,Q,X0,P0,R,xyz,Z)     
    
    Xn=F*X0;%预测
    vec=[Xn(1);Xn(4);Xn(7)]-xyz;
    
    %象限获取
    if vec(1)>=0
        if vec(2)>=0
            dem=0;%第一、四象限内不需要对角度做修正
        else
            dem=0;
        end       
    else
         if vec(2)>=0%第二象限要加pi
            dem=1;
         else
            dem=-1;%第三象限要减pi
         end  
    end
    
    Zt=[sqrt((vec(1))^2+(vec(2))^2+(vec(3))^2);...
        atan((vec(2))/(vec(1)));...
        asin((vec(3))/(sqrt((vec(1))^2+(vec(2))^2+(vec(3))^2)))];
    
    Zt(2)=Zt(2)+pi*dem;%获得自由量测
    
    dh1_dx= vec(1)/sqrt(vec(1)^2+vec(2)^2+vec(3)^2);
    dh1_dy= vec(2)/sqrt(vec(1)^2+vec(2)^2+vec(3)^2);
    dh1_dz= vec(3)/sqrt(vec(1)^2+vec(2)^2+vec(3)^2);
    
    dh2_dx=-vec(2)/(vec(1)^2+vec(2)^2);
    dh2_dy= vec(1)/(vec(1)^2+vec(2)^2);
    dh2_dz= 0;
    
    dh3_dx=-vec(1)*vec(3)/(vec(1)^2+vec(2)^2+vec(3)^2)/sqrt(vec(1)^2+vec(2)^2);
    dh3_dy=-vec(2)*vec(3)/(vec(1)^2+vec(2)^2+vec(3)^2)/sqrt(vec(1)^2+vec(2)^2);
    dh3_dz=sqrt(vec(1)^2+vec(2)^2)/(vec(1)^2+vec(2)^2+vec(3)^2);
%     dh3_dx=-vec(1)*vec(3)/sqrt(vec(1)^2+vec(2)^2+vec(3)^2)/(vec(1)^2+vec(2)^2);    
%     dh3_dy=-vec(2)*vec(3)/sqrt(vec(1)^2+vec(2)^2+vec(3)^2)/(vec(1)^2+vec(2)^2);
        
    H=[dh1_dx,0,0,dh1_dy,0,0,dh1_dz,0,0;
       dh2_dx,0,0,dh2_dy,0,0,dh2_dz,0,0;
       dh3_dx,0,0,dh3_dy,0,0,dh3_dz,0,0];
  
    P=F*P0*F'+Q;
    Pzminus=H*P*H'+R;
    K=P*H'/Pzminus;
    v=Z-Zt;
    X0=Xn+K*(Z-Zt);
    P0=(eye(9)-K*H)*P;
end

