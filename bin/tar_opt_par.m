function [ R,mod_dif ] = tar_opt_par( p,tg_r,res_sig,wininfo,psize,we,ntol)
%para to matrix form
P=reshape(p,psize(1),psize(2));

%evaluate model on wide time grid
mod_neg=zeros(size(tg_r));
for i=1:psize(2)
    mod_neg=mod_neg+mod_ga(tg_r,P(:,i));
end
%extract model on narrow grid
mod_dif=mod_neg(wininfo(1):wininfo(1)+wininfo(2)-1);

%get residual
mod_resi=res_sig(wininfo(1):wininfo(1)+wininfo(2)-1)-mod_dif;

if ntol>0
    id=abs(mod_resi)<=ntol;            
            
    mod_resi(id)=0;
end

% model fit on narrow window
R1=we(1)*(mod_resi);

% model nonnegativity on wide window
R2=we(2)*0.1*min(res_sig-mod_neg,0);

%combine residuals
R=[R1(:); R2(:)];
end

