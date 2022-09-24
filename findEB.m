function [fitness] = findEB(KRij,Rij,Req)
Sum=0;    
b=length(KRij);
for i=1:b
    Sum=Sum+(KRij(i)/2)*(Rij(i)-Req(i))^2;
end
fitness=Sum;
end