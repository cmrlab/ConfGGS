function t_deg = findAngle(a,b,c)
    ba=a-b;
    bc=c-b;
    modba=findmod(ba);
    modbc=findmod(bc);
    n=sum(ba.*bc);
    t_rad=acos(n/(modba*modbc));
    t_deg=(t_rad*180)/3.141;
end