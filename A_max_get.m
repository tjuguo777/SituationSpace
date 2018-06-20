function [a_max] = A_max_get(ngz,thea)
    g=9.8;
    ax=g*(ngz(2)^2-1)^0.5;
    ay=g*(ngz(1)-sin(thea));
    az=g*(ngz(2)-cos(thea));
    a_max=[ax;ay;az];
end

