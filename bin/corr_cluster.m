function [rescl] = corr_cluster(gx,gc,para)
C_len=para.cl_n;
NoiseWin=para.NoiseWin;
gcs=gc;


test=zeros(1,length(gcs)-NoiseWin+1);
for i=1:length(test)
    gctest=gc(i:i+NoiseWin-1);
    switch 'maxmin'
        case 'maxmin'
            test(i)=max(gctest)-min(gctest);
        case 'std'
            test(i)=6*std(gctest);
    end
end
[mi,idx]=min(test);


% set noise level
C_tol=para.cl_fac*mi;


C_XGidx=round(linspace(1,length(gx),C_len));
C_gcs=cell(C_len,1);
C_gx=cell(C_len,1);
C_act=false(C_len-1,1);
C_act_x=zeros(C_len-1,1);
for i=1:C_len-1
    C_D=gcs(C_XGidx(i):C_XGidx(i+1)-1);
    C_gcs{i}=C_D;
    C_gx{i}=gx(C_XGidx(i):C_XGidx(i+1)-1);
    C_temp=linspace(0,1,length(C_D));
    A=[ones(length(C_temp),1)];
    aa=A\C_D;
    C_Dc=C_D-A*aa;
    if min(abs(C_Dc)<=C_tol)==1
        C_act(i)=1;
        C_act_x(i)=gx(round(1/2*(C_XGidx(i)+C_XGidx(i+1)-1)));
    end
end

C_act=logical(floor(smooth(C_act+0.0,5)));
rescl.C_act_x=C_act_x;
rescl.C_act=C_act;

% get sub systems
act=~C_act;
targets=[];
area_ct=0;
for i=1:length(act)
    if i==1
        if act(i)==1
            area_ct=1;
            t1=1;
        end
    end
    if i>1
        if act(i-1)==0 && act(i)==1
            area_ct=area_ct+1;
            t1=i;
        end
        if act(i-1)==1 && act(i)==0
            targets=[targets;[t1 i]];
        end
    end
    if i==length(act)
        if act(i)==1
            targets=[targets;[t1 i]];
        end
    end
end

Cgc=cell(area_ct,1);
Cgx=cell(area_ct,1);
for i=1:area_ct
    t1=C_XGidx(targets(i,1));
    t2=C_XGidx(targets(i,2));
    
    Cgc{i}=gc(t1:t2);
    Cgx{i}=gx(t1:t2);
end


rescl.anz_prob=area_ct;
rescl.Cgc=Cgc;
rescl.Cgx=Cgx;
rescl.noiselvl=C_tol;






end

