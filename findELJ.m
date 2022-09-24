function [fitness] = findELJ(Aij,Bij,Rij)
Sum=0;    
n=length(Aij);
    for i=1:n
        Sum=Sum+(Aij(i)/Rij(i))^12-(Bij(i)/Rij(i))^6;
    end
fitness=Sum;
end