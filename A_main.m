%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code modified:2018/5/29                                                  %
%function     :real-time path planning of UCAV                            %
%athours      :ysx---youshixun@hrbeu.edu.cn                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        A. objects initialization                        %
%                                                                         %
%-------------------------------------------------------------------------%

Lo=[30,60,90,450];                               %������뾶Ϊ30 60 90 330 m
Lu=30;                                                       %�ɻ��뾶Ϊ30 m

global v_max Num_objects T C a_max
v_max=250;                                              % UCAV����ٶ�250m/s
T=0.5;                                                            % ��������

NT=200;                                                             % ������
Model=3;                                                       % IMMģ������
ngz=[8,8];                                           % ���������غͷ������
muda_att=[2 4 8 12];                                              % ��Ӧϵ��

%----------------Optional parameters and environment----------------------%
muda=muda_att(4);                                             % ѡ����Ӧϵ��
Scenario=1;                                                   % ѡ����Ի���
K_show_for_threats=[40 60;85 105;130 150];                        % ��ʾλ��
K=200;                                                          % ʱ������
filter_mode=2;                                                 % �˲�ģʽCKF
%-----------------------------for testing---------------------------------%

UCAV(:,1)=[15e3;30;0;2e3;0;0;6e3;0;0];                      % UCAV����ʼ״̬
UCAV_p(:,1)=[UCAV(1,1);UCAV(4,1);UCAV(7,1)];                      % ��ʼλ��
UCAV_v(:,1)=[UCAV(2,1);UCAV(5,1);UCAV(8,1)];                      % ��ʼ�ٶ�
UCAV_a(:,1)=[UCAV(3,1);UCAV(6,1);UCAV(9,1)];                    % ��ʼ���ٶ�

if Scenario==1
    X{1,1}=[2e3;80;0;2e3;-30;0;5e3;-60;0];                        % ����Ŀ��
else
    X{1,1}=[2e3;60;0;1e3;90;0;4e3;50;0];               % 1������Ŀ��3���ϰ���
    X{2,1}=[6e3;57;6.5;4e3;-65;0;2.5e3;140;-3]; 
    X{3,1}=[1.6e3;90;0;4e3;0;0;7.8e3;-15;0];
    X{4,1}=[5e3;0;0;3e3;0;0;5.5e3;0;0];                          % ��ֹ����
end

Num_objects=size(X,1);                                              % Ŀ����
p=[4e4,0,0;                                                     % ��ʼЭ����
    0,16,0;
    0,0,1];
w=pi*[9/180 3/180 -9/180 0];                                  % ��ת�Ƕ���Ϣ

for n_o=1:Num_objects                                             % �˶�����
    [F{1,n_o},Q{1,n_o}]=Motion_cvf(T);
    [F{2,n_o},Q{2,n_o}]=Motion_caf(T);
    [F{3,n_o},Q{3,n_o}]=Motion_ctf(w(n_o),T); 
    if Scenario==1
        [F{3,n_o},Q{3,n_o}]=Motion_ctf(w(2),T);
    end   
end

for n_o=1:Num_objects                                           % ��ʼ��״̬
    X_immckf{n_o,1}=X{n_o,1};                                 
    core_x{n_o,1}=[X_immckf{n_o,1}(1);
                   X_immckf{n_o,1}(4);
                   X_immckf{n_o,1}(7)];                 
    vx{n_o,1}=[X_immckf{n_o,1}(2);X_immckf{n_o,1}(5);X_immckf{n_o,1}(8)];
    ax{n_o,1}=[X_immckf{n_o,1}(3);X_immckf{n_o,1}(6);X_immckf{n_o,1}(9)];   
    P{n_o,1}=blkdiag(p,p,p);                  
    u{n_o,1}=[0.8,0.1,0.1];           % ��ʼ�˶����ʷֲ�[���� �ȼ��� ƽ����ת]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           B. maneuvering                                %
%                                                                         %
%-------------------------------------------------------------------------%

for k=1:NT
    
%-------------------------------------------------------------------------%
%                          1. equation of motion                          %
%-------------------------------------------------------------------------%

    if Scenario==1
        X{1,k+1}=F{3,1}*X{1,k}+sqrtm(Q{3,1})*randn(9,1); 
    else
        if k<31
            X{1,k+1}=F{1,1}*X{1,k}+sqrtm(Q{1,1})*randn(9,1);      % ����Ŀ��
        else
            if k<101
                X{1,k+1}=F{3,2}*X{1,k}+sqrtm(Q{3,2})*randn(9,1);
            else 
                X{1,k+1}=F{3,1}*X{1,k}+sqrtm(Q{3,1})*randn(9,1);
            end
        end
        if k<81
            X{2,k+1}=F{2,2}*X{2,k}+sqrtm(Q{2,2})*randn(9,1);       % �ϰ���1
        else
            X{2,k+1}=X{2,k};
        end
        X{3,k+1}=F{3,3}*X{3,k}+sqrtm(Q{3,3})*randn(9,1);           % �ϰ���2
        X{4,k+1}=F{1,4}*X{4,k}+sqrtm(Q{1,4})*randn(9,1);           % �ϰ���3
    end
    
%-------------------------------------------------------------------------%
%                         2. measurement generation                       %
%-------------------------------------------------------------------------%
    
    target(:,k)=[X{1,k}(1);X{1,k}(4);X{1,k}(7)];
    
    if  Scenario==2
        threat1(:,k)=[X{2,k}(1);X{2,k}(4);X{2,k}(7)];
        threat2(:,k)=[X{3,k}(1);X{3,k}(4);X{3,k}(7)];
        threat3(:,k)=[X{4,k}(1);X{4,k}(4);X{4,k}(7)];
    end
    
    for n_o=1:Num_objects               % ����������ݣ��Լ�����������������
        [Z_x{n_o,k},oxyz{n_o,k},R{n_o,k}]=Z_production(X{n_o,k},UCAV_p(:,k));
    end
    
%-------------------------------------------------------------------------%
%                                3. filtering                             %
%-------------------------------------------------------------------------%
    
    if k==1
        for n_o=1:Num_objects
            for i=1:Model
                Xn_est{i,n_o}=X_immckf{n_o,1};
                Pn_est{i,n_o}=P{n_o,1};
            end
            Xn_nl{n_o,1}=IMM_prediction(F,Q,Xn_est,...
                                        Pn_est,u{n_o,1},n_o);       % Ԥ��ֵ
        end
    else
        for n_o=1:Num_objects                                          
            [X_immckf{n_o,k},P{n_o,k},Xn_est,Pn_est,u{n_o,k}]=...
            IMM_filtering(F,Q,Z_x{n_o,k},R{n_o,k},UCAV_p(:,k),...
                          Xn_est,Pn_est,u{n_o,k-1},n_o,filter_mode);            % �˲�ֵ
                      
            Xn_nl{n_o,k}=IMM_prediction(F,Q,Xn_est,...
                                        Pn_est,u{n_o,k},n_o);       % Ԥ��ֵ
                                    
            core_x{n_o,k}=[X_immckf{n_o,k}(1);
                           X_immckf{n_o,k}(4);
                           X_immckf{n_o,k}(7)];                   % �˲�λ��
                                    
            if n_o==3                                      
                vx{n_o,k}=[-X_immckf{n_o,k}(2);
                           -X_immckf{n_o,k}(5);
                            X_immckf{n_o,k}(8)];                % �˲����ٶ�
            else
                vx{n_o,k}=[X_immckf{n_o,k}(2);
                           X_immckf{n_o,k}(5);
                           X_immckf{n_o,k}(8)];
            end
            
            ax{n_o,k}=[X_immckf{n_o,k}(3);
                       X_immckf{n_o,k}(6);
                       X_immckf{n_o,k}(9)];                   % �˲��ļ��ٶ�
                                                                  
            dr(n_o,k)=norm([X{n_o,k}(1);X{n_o,k}(4);X{n_o,k}(7)]...
                           -UCAV_p(:,k));                         % ��ʵ����
        end
    end
    
    [C,a_max]=XYZ_frame_change(UCAV_v(:,k),ngz);    % ��ء��ɻ�������ϵת��,
                                                        % �Լ���Ӧ��������
     
%-------------------------------------------------------------------------%
%                            4. situation space                           %
%-------------------------------------------------------------------------%
    
    UCAV_nl=F{2,1}*UCAV(:,k);                               % UCAV��Ԥ��״̬
    core_u=UCAV_p(:,k);                                          % UCAVԲ��
    ru=Lu+norm([UCAV_nl(1);UCAV_nl(4);UCAV_nl(7)]-core_u);         % �뾶ru
    
    for n_o=1:Num_objects
        rt(n_o)=Lo(n_o)+norm([Xn_nl{n_o,k}(1);
                              Xn_nl{n_o,k}(4);
                              Xn_nl{n_o,k}(7)]...
                             -core_x{n_o,k});               % Ŀ��İ�ȫ�뾶 
                                                                        
        dvt(n_o)=norm([X_immckf{n_o,k}(1);
                       X_immckf{n_o,k}(4);
                       X_immckf{n_o,k}(7)]-UCAV_p(:,k));      
                   
        LOS(:,n_o)=([X_immckf{n_o,k}(1);
                     X_immckf{n_o,k}(4);
                     X_immckf{n_o,k}(7)]...
                    -UCAV_p(:,k))/dvt(n_o);                   % Ŀ�����߷���
        
        vut(:,n_o)=UCAV_v(:,k)-vx{n_o,k}; 
        vut_p=vut(:,n_o)'*LOS(:,n_o);                         % ͶӰ����ٶ�
        a_p=LOS'*ax{1,k};                                   % ͶӰ��Լ��ٶ�
        a_gen=norm(Constraint_boundary(LOS(:,n_o),a_max));  % ���ۼ��ٶ��Ͻ�
        
        rou=vut_p^2/2/norm(-a_gen-a_p);                       % ��ȫ�������
        rx(n_o,k)=rt(n_o)+ru+rou;                                 % ̬�ư뾶
    end
    
%-------------------------------------------------------------------------%
%                            5. target tracking                           %
%-------------------------------------------------------------------------%  

    at_expcet=Attractive(dvt(1),rx(1,k),...
                         LOS(:,1),vut(:,1),muda);               % ���ټ��ٶ�
    
%-------------------------------------------------------------------------%
%                           6. collision avoidance                        %
%-------------------------------------------------------------------------%
    
    ab_expcet=Collision(dvt,rx(:,k),LOS,vut,muda);              % ���ϼ��ٶ�

%-------------------------------------------------------------------------%
%                             7. path planning                            %
%-------------------------------------------------------------------------%

    if norm(ab_expcet)==0
        aplan=Repair(at_expcet,UCAV_v(:,k));                % ���ٶȽ����޸�
    else
        aplan=Repair(ab_expcet,UCAV_v(:,k));
    end
    
    UCAV(3,k)=aplan(1);UCAV(6,k)=aplan(2);UCAV(9,k)=aplan(3);
    UCAV(:,k+1)=F{2,1}*UCAV(:,k);
    UCAV_p(:,k+1)=[UCAV(1,k+1);UCAV(4,k+1);UCAV(7,k+1)];
    UCAV_v(:,k+1)=[UCAV(2,k+1);UCAV(5,k+1);UCAV(8,k+1)];      % �����µ�״̬

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        C. path visualization                            %
%                                                                         %
%-------------------------------------------------------------------------%

% figure
% hold on;
% grid on;
set(gca,'gridlinestyle','--');
% set(gcf,'color','white');
% set(gca,'tickdir','in')
% 
% Draw_sphere(3*Lu/1e3,UCAV_p(:,1)/1e3, [0,0,0]);            % ��ɫ-UCAV��ʼ��
% xlabel('x/km');ylabel('y/km');zlabel('z/km')
% 
% for k=1:K
%     Draw_sphere(Lu/1e3,UCAV_p(:,k)/1e3, [0,0,1]);               % ��ɫ-UCAV
%     Draw_sphere(Lo(1)/1e3,target(:,k)/1e3, [1,0,0]);             % ��ɫ-Ŀ��
%     if Scenario==2
%         axis([0,16,0,7,1,10]);
% %         if k>=K_show_for_threats(3,1)&&k<=K_show_for_threats(3,2)
%             Draw_sphere(Lo(2)/1e3,threat1(:,k)/1e3,...
%                         [1,0,1]);                           % õ��ɫ-���ϰ�1
%             Draw_sphere(Lo(3)/1e3,threat2(:,k)/1e3,... 
%                         [0.2,0.8,0]);                         % ��ɫ-���ϰ�2
%             Draw_sphere(Lo(4)/1e3,threat3(:,k)/1e3,...
%                         [0.5,0.5,0.5]);                      % ��ɫ-��̬�ϰ�
% %         end
%     else
%         axis([0,16,0,6,0,6]);
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          D. path reliability                            %
%                                                                         %
%-------------------------------------------------------------------------%

% figure
subplot(3,1,1)
box on;
hold on;
plot(1:K,zeros(1,K),'k--');
set(gcf,'color','white');
plot(1:K,log(dr(1,:)/(Lu+Lo(1))),'r');

if Scenario==2
    plot(1:K,log(dr(2,:)/(Lu+Lo(2))),'m');
    plot(1:K,log(dr(3,:)/(Lu+Lo(3))),'c');
    plot(1:K,log(dr(4,:)/(Lu+Lo(4))),'g');
end

axis([1 NT 0 6])
xlabel('Sampling point');ylabel('Absolute safety factor')      %���԰�ȫϵ��

subplot(3,1,2)
box on;
plot(1:K,zeros(1,K),'k--');
hold on;
plot(1:K,log(dr(1,:)./rx(1,:)),'r');

if Scenario==2
    plot(1:K,log(dr(2,:)./rx(2,:)),'m');
    plot(1:K,log(dr(3,:)./rx(3,:)),'c');
    plot(1:K,log(dr(4,:)./rx(4,:)),'g');
end

axis([1 NT -1 5])
xlabel('Sampling point');ylabel('Relative safety factor')      %��԰�ȫϵ��

subplot(3,1,3)
box on;
hold on;
plot(1:K,sum(abs(UCAV_v(:,1:K)).^2).^0.5,'m');

if Scenario==1
    plot(1:K,104*ones(1,K),'k--');
else
    plot(1:K,119*ones(1,K),'k--');
end

axis([1 NT 0 300])
set(gcf,'color','white');  
xlabel('Sampling point');ylabel('Speed (m/s)')                 %�ٶ���������