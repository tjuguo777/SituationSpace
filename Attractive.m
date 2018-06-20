function [act_expcet]=Attractive(dvt,rx,LOS,vut,bl)
    global T C a_max
    
%-------------------------------------------------------------------------%
%                           1. adjustment of v                            %
%-------------------------------------------------------------------------%

    vtu=-vut;
    avc=C*(vtu/T);                    
    dav=Constraint_boundary(avc/norm(avc),a_max);   
    if norm(avc)>norm(dav)
         avc=dav;
    end
    
%-------------------------------------------------------------------------%
%                          2. adjustment of dis                           %
%-------------------------------------------------------------------------%

    as1=2*(dvt-rx-vut'*LOS*T)/(T^2)*LOS;
    as2=-(vut-vut'*LOS*LOS)/T;
    as=as1+as2;             
    asc=C*as;                        
    as_cb=Constraint_boundary(asc/norm(asc),a_max);
    if norm(asc)>norm(as_cb)
         asc=as_cb;
    end
    
%-------------------------------------------------------------------------%
%                          3. adaptation synthesis                        %
%-------------------------------------------------------------------------%
    
    if dvt<bl*rx&&dvt>rx
        cs=(dvt/rx-1)/(bl-1);
    else
        cs=1;
    end
    cv=1-cs;
    act_expcet=cv*avc+cs*asc; 
    if norm(Constraint_boundary(act_expcet/norm(act_expcet),a_max))...
                                                  <norm(act_expcet)
        act_expcet=Constraint_boundary(act_expcet/norm(act_expcet),a_max);
    end
    
%-------------------------------------------------------------------------%
%                       4. projection transformation                      %
%-------------------------------------------------------------------------%

    act_expcet=C'*act_expcet;
end