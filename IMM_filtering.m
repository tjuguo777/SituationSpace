function [X_est,P_est,Xn_est,Pn_est,u]=IMM_filtering(F,Q,Z,R,xyz0,...
                                                     Xn_est,Pn_est,u,n_o,mode)

%-------------------------------------------------------------------------%
%                          1. IMM initialization                          %
%-------------------------------------------------------------------------%
    model=size(u,2);
    Pr=[0.8,0.1,0.1;
        0.1,0.8,0.1;
        0.1,0.1,0.8];                                          %转移概率矩阵 
    c_mean=zeros(1,model);                                       %归一化常数
    mu=zeros(model,model);
    
    for i=1:model
        c_mean=c_mean+Pr(i,:)*u(i);
    end
    
    for i=1:model
        mu(i,:)=Pr(i,:)*u(i)./c_mean;                        %贝叶斯转移概率
    end
    
    for j=1:model
        X0{j,1}=zeros(9,1);
        P0{j,1}=zeros(9);
        for i=1:model
            X0{j,1}=X0{j,1}+Xn_est{i,n_o}*mu(i,j);             %混合状态估计
        end
        for i=1:model
            P0{j,1}=P0{j,1}+mu(i,j)*( Pn_est{i,n_o}...       %混合协方差估计
                    +(Xn_est{i,n_o}-X0{j,1})*(Xn_est{i,n_o}-X0{j,1})');
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 这部分用了一个小trick，因为MATLAB的数据精度有限，滤波在迭代时可能出现精度过小%
% 当似然度精度小于1e-320时，分母将全部置0，故提前对贝叶斯概率进行分子相消，    %
% 保证分母有常数                                                           %                                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    % A=zeros(1,model);%似然矩阵
    A_matrix=eye(model);%似然矩阵
   
    for j=1:model 
                        
%-------------------------------------------------------------------------%
%                             2. IMM+EKF/CKF/UKF                          %
%-------------------------------------------------------------------------% 
    if mode==1        
        [Xn_est{j,n_o},Pn_est{j,n_o},v{j,1},s{j,1}]=...      %各自通道的滤波
                            EKF_filtering_model(F{j,n_o},Q{j,n_o},...
                                                X0{j,1},P0{j,1},R,xyz0,Z);                        
    end

    if mode==2
        [Xn_n{j,1},Pminus]=CKF_prediction_model(F{j,n_o},Q{j,n_o},...
                                                X0{j,1},P0{j,1});     %预测
                                                                                        
        [Xn_est{j,n_o},Pn_est{j,n_o},v{j,1},s{j,1}]=...      %各自通道的滤波
                            CKF_filtering_model(Xn_n{j,1},Pminus,R,xyz0,Z); 
    end
    
    if mode==3
        [Xn_est{j,n_o},Pn_est{j,n_o},v{j,1},s{j,1}]=...      %各自通道的滤波
                            UKF_filtering_model(F{j,n_o},Q{j,n_o},...
                                                X0{j,1},P0{j,1},R,xyz0,Z);     
    end
    
    
%-------------------------------------------------------------------------%
%                           2. update the output                          %
%-------------------------------------------------------------------------%  

        % nl=length(s{j,1})/2;
        % A(j)=1/((2*pi)^nl*sqrt(det(s{j,1})))*exp(-0.5*v{j,1}'/(s{j,1})*v{j,1});
        for i=1:j
            if i~=j           
                A_matrix(i,j)=sqrt(det(s{i,1})/det(s{j,1}))*...
                                   exp(0.5*(v{i,1}'/...
                                 (s{i,1})*v{i,1}-v{j,1}'/(s{j,1})*v{j,1}));
                A_matrix(j,i)=sqrt(det(s{j,1})/det(s{i,1}))*...
                                   exp(0.5*(v{j,1}'/...
                                 (s{j,1})*v{j,1}-v{i,1}'/(s{i,1})*v{i,1}));
            end
        end      
    end  
    
    u=c_mean./(c_mean*A_matrix'); 
%     u=(A.*c_mean)/sum(A.*c_mean);                                
   
    Xn=zeros(9,1);
    Pn=zeros(9);

    for j=1:model
        Xn=Xn+Xn_est{j,n_o}.*u(j);
    end
    
    for j=1:model
        Pn=Pn+u(j).*(Pn_est{j,n_o}+(Xn_est{j,n_o}-Xn)*(Xn_est{j,n_o}-Xn)');
    end
    
    X_est=Xn;                                                  %统计状态输出
    P_est=Pn;                                                %统计协方差输出
end

