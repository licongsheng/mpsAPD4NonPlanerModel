function in_tri = is_in_tri(P,v_tri)
    A = v_tri(1,:);
    B = v_tri(2,:);
    C = v_tri(3,:);
    
    AP = P-A;
    x1 = AP(1);
    y1 = AP(2);
    AC = C-A;
    x2 = AC(1);
    y2 = AC(2);
    AB = B-A;
    x3 = AB(1);
    y3 = AB(2);
    
    m = x2*y3-x3*y2;
    u = (x1*y3-x3*y1)/m;
    v = -(x1*y2-x2*y1)/m;
    
    if(u>=0&&v>=0&&u+v<=1)
        in_tri = [u,v];
    else
        in_tri = 0;
    end
end