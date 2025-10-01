%calculate the Frobenius form of c1
function [comp]=c1f(ifim)
[u,s,v]=svd(ifim);
p=rank(ifim);
s=diag(s);
s=s(1:p);
sbar=sum(s)/p;
comp=(s-sbar)'*(s-sbar)/(4*sbar^2);
return