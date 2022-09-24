function s=check(T,x,y)
s=1;
    for i=1:size(T,1)
        if(T(i,1)==x && T(i,2)==y)
            s=0;
            break;
        end
    end
end