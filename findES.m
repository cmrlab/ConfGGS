function [fitness] = findES(qi,qj,Rij,chro)
e=1;
% Sum=qi*qj/(e*Rij);
% fitness=Sum;
sum=0;
charge=csvread('charge.csv');
cod_bond=csvread('cod_bond.csv');
[x y]=size(cod_bond);
T=zeros(x*2,2);
T(1:x,1:2)=cod_bond(:,1:2);
T(x+1:2*x,1)=cod_bond(:,2);
T(x+1:2*x,2)=cod_bond(:,1);
T=sortrows(T,1);
for i=1:7
    for j=i+1:8
        if check(T,i,j)==1
            sum=sum+(charge(i)*charge(j))/(e*(findDist(chro(i,:),chro(j,:))));
        end
    end
end
fitness=sum;
end