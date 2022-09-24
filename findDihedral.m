function t_deg = findDihedral(p1,p2,p3,p4)
    q1=p2-p1;
    q2=p3-p2;
    q3=p4-p3;

    cq1q2=cross(q1,q2);
    cq2q3=cross(q2,q3);
    n1=cross(q1,q2)/absolute(cross(q1,q2));
    n2=cross(q2,q3)/absolute(cross(q2,q3));
    u1=n2;
    u3=q2/absolute(q2);
    u2=cross(u3,u1);
    cost=dot(n1,u1);
    sint=dot(n1,u2);
    t_rad=-atan2(sint,cost);
    t_deg=(t_rad*180)/3.141;
end

