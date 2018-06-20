function [ Fct,Qct] = Motion_ctf(w,T)
    w=w+1e-6;
    Qcts=[(6*w*T-8*sin(w*T)+sin(2*w*T))/4/w^5,2*sin(0.5*w*T)^4/w^4,(-2*w*T+4*sin(w*T)-sin(2*w*T))/4/w^3;
        2*sin(0.5*w*T)^4/w^4,(2*w*T-sin(2*w*T))/4/w^3,sin(w*T)^2/2/w^2;
       (-2*w*T+4*sin(w*T)-sin(2*w*T))/4/w^3,sin(w*T)^2/2/w^2,(2*w*T+sin(2*w*T))/4/w];
    Qct=1e-4*blkdiag(Qcts,Qcts,Qcts);%×´Ì¬Îó²îÐ­·½²î
    Fcts=[1 sin(w*T)/abs(w) 0 0 -(1-cos(w*T))/abs(w) 0;
          0 cos(w*T)        0 0 -sin(w*T)            0;
          0 -w*sin(w*T)     0 0 -w*cos(w*T)          0;
          0 (1-cos(w*T))/abs(w) 0 1 sin(w*T)/abs(w)  0;
          0 sin(w*T)        0 0      cos(w*T)        0;
          0  w*cos(w*T)     0 0       -w*sin(w*T)    0];
    Fcvs=[1 T 0;
          0 1 0
          0 0 0];
    Fct=blkdiag(Fcts,Fcvs);%×´Ì¬×ªÒÆ¾ØÕó
end

