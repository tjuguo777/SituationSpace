function [X0,P0,v,Pzz]=UKF_filtering_model(F,Q,X0,P0,R,xyz,Z)

    L=numel(X0);                                            % UT�任ϵ������
    m=numel(Z);
    alpha=1;
    ki=0;
    beta=2.5;
    rambda=alpha^2*(L+ki)-L;
    c=L+rambda;
    
    Wm=[rambda/c 0.5/c+zeros(1,2*L)];                   % UT�任��Ȩϵ������
    Wc=Wm;
    Wc(1)=Wc(1)+(1-alpha^2+beta);                % ��Э�����������Э���ͬ
    c=sqrt(c);
    Xsigmaset=sigmas(X0,P0,c);                              % UT�任�Ĺ����
    l=size(Xsigmaset,2);

    dem=zeros(1,l);
    f=@(x)(F*x);
    [Xmeans,Xsigma_pre,Pn,Xdiv]=ut(f,Xsigmaset,dem,...
                                   Wm,Wc,L,Q);%����Ԥ��ֵ����ֵ��״̬��Э����
      
    for i=1:l
        vec(1:3,i)=[Xsigma_pre(1,i);
                    Xsigma_pre(4,i);
                    Xsigma_pre(7,i)]-xyz;
    
        if vec(1,i)>=0
                dem(i)=0;    
        else
            if vec(2,i)>=0
               dem(i)=1;
            else
               dem(i)=-1;
            end  
        end
    end
  
    h=@(x)[sqrt((x(1))^2+(x(2))^2+(x(3))^2);...
           atan((x(2))/(x(1)));...
           asin((x(3))/(sqrt((x(1))^2+(x(2))^2+(x(3))^2)))];
   
    [Zpre,Z_,Pzz,Zdiv]=ut(h,vec,dem,Wm,Wc,m,R);
    
    Pxz=Xdiv*diag(Wc)*Zdiv';
    K=Pxz/Pzz;
    v=Z-Zpre;
    X0=Xmeans+K*v;
    P0=Pn-K*Pxz';   
end

