function [fitness] = findEA(KTijk,Tijk,Teq)
Sum=0;    
a=length(KTijk);
for i=1:a
    Sum=Sum+(KTijk(i)/2)*(Tijk(i)-Teq(i))^2;
end
fitness=Sum;
end

