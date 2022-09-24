function [fitness] = findED(Vmijkl,Fiijkl,Gijkl)
oSum=0;    
m=length(Vmijkl);
d=m;
for j=1:d
    inSum=0;
    for i=1:m
        inSum=inSum+(Vmijkl(i)/2)*(1+cos(i*Fiijkl(i)-Gijkl(i)));
    end
    oSum=oSum+inSum;
end
fitness=oSum;
end

