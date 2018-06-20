function h = Draw_sphere(r, center, color)
    [x,y,z] = sphere(50);
    h = surf(r*x+center(1), r*y+center(2), r*z+center(3));
    if isempty(color)
        h.EdgeColor = rand(1,3);
        h.FaceColor = h.EdgeColor;
    else
        h.EdgeColor = color;
        h.FaceColor = color;   
    end
    axis equal
end