function [ Fcv,Qcv] =Motion_cvf(T)
    Fcvs=[1 T 0;
          0 1 0
          0 0 0];
    Fcv=blkdiag(Fcvs,Fcvs,Fcvs);%×´Ì¬×ªÒÆ¾ØÕó
    Qcvs=[1,0,0;
          0,0.01,0;
          0,0,1e-4];
    Qcv=blkdiag(Qcvs,Qcvs,Qcvs);%×´Ì¬Îó²îÐ­·½²î
end

