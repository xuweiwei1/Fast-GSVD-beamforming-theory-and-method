function [D1,D2,W]=fgsvdbeamforming(H1,H2)
[U,V,X]=rgsvd4(H1,H2);
D1=U';
D2=V';
[u,v,c]=svd(X);
W=c*pinv(v)*u';
end