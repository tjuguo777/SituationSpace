function [ Fcv,Qcv] =Motion_cvf(T)
    Fcvs=[1 T 0;
          0 1 0
          0 0 0];
    Fcv=blkdiag(Fcvs,Fcvs,Fcvs);%״̬ת�ƾ���
    Qcvs=[1,0,0;
          0,0.01,0;
          0,0,1e-4];
    Qcv=blkdiag(Qcvs,Qcvs,Qcvs);%״̬���Э����
end

