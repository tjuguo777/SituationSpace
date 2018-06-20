function [ Fca,Qca] =Motion_caf(T)
    Fcas=[1 T T^2/2;
          0 1 T
          0 0 1];
    Fca=blkdiag(Fcas,Fcas,Fcas);%×´Ì¬×ªÒÆ¾ØÕó
    Qcas=[1,0,0;
          0,0.01,0;
          0,0,1e-4];
    Qca=blkdiag(Qcas,Qcas,Qcas);%×´Ì¬Îó²îÐ­·½²î
end

