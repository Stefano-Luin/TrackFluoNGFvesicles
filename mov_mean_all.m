function out=mov_mean_all(v,tc,dt)
% Copyright (C) 2012 - 2022 Carmine di Rienzo; Stefano Luin (s.luin@sns.it)
l=length(v);
if nargin<3 || isempty(dt)
    dt=ones(size(v));
end
if nargin<2 || isempty(tc) || tc<=1
    out=v;
    return;
end
tc=ceil(tc);
if l<tc
    out=(v(:)'*dt(:))/sum(dt(:));
else
    vv=v(:);dtt=dt(:);
    out=arrayfun(@(c) (vv(c:(c+tc-1))'*dtt(c:(c+tc-1)))/sum(dtt(c:(c+tc-1))),1:l-tc+1);
end