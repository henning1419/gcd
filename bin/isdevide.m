function [b] = isdevide(gr,sig,minlen,relmaxpeakW)
countermax=80;

if length(gr)<minlen
    b=0;
    return
end
len=gr(end)-gr(1);

sig=smooth(sig,40,'sgolay',3);


[ma,id]=max(sig);
id=max(id,2);
id=min(id,length(gr)-1);

le=id;
ri=id;

cc=0;
while cc<countermax
    if le==1
        break
    end
    while (sig(le)>=sig(le-1) && le>1)
        le=le-1;
        if le==1
            break
        end
    end
    if le>1
        le=le-1;
    end
    cc=cc+1;
end

cc=0;
while cc<countermax
    if ri==length(gr)
        break
    end
    while sig(ri)>=sig(ri+1) && ri<=length(gr)-1
        ri=ri+1;
        if ri==length(gr)
            break
        end
    end
    if ri<=length(gr)-1
        ri=ri+1;
    end
    cc=cc+1;
end

if (gr(ri)-gr(le))<relmaxpeakW*len
    b=1;
else
    b=0;
end

end

