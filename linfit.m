function [p,cov,ndf,Chi2red,Chi2]=linfit(x,y,w)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
if nargin==2
    [p,S]=polyfit(x,y,1);
    if nargout>1, cov=covpolyfit(S); end
else
    if nargin<3 || isempty(w)
        w=ones(size(x));
    end
    x=reshape(x,[],1);
    y=reshape(y,[],1);
    w=sqrt(reshape(w,[],1));
    X=[w,w.*x];
    Y=w.*y;
    Q=X'*X;
%     M=inv(Q);
    p=Q\(X'*Y);
%     p=M*(X'*Y);
    if nargout>1
        ndf=(length(x)-2);
        Chi2=sum((Y-X*p).^2);
        Chi2red=Chi2/ndf;
        cov=Chi2red*inv(Q);
%          cov=Chi2red*M;
    end
end
return

