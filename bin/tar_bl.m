function [ R,mod ] = tar_bl( hi,x,y,center)

mod = spline(center,hi,x);
% mod = interp1(center,hi,x,'linear');

R1=10*min(y-mod,0);
R2=0.0001*abs(1-mod);

R=[R1(:);R2(:)];

end

