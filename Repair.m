function [aplan]=Repair(a,vu)
    global v_max T
        k=norm(a*T+vu)/v_max;
    if k>1
        aplan=((a*T+vu)/k-vu)/T;
    else
        aplan=a;
    end
end