function x=absolute(y)
sum=0;
for i=1:length(y)
    sum=sum+y(i)^2;
end
x=sqrt(sum);
end