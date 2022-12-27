function [byt,celldata]=scrivi(fname,data,txt,header,sep)
%byt=SCRIVI(fname,data,txt,header,sep)
%   see also APRI WRITEFILE READFILE OPENASCII
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
if nargin<5, sep=''; end
sizes=size(data);
if nargin<3 || isempty(txt)
    txt={};
    add=0;
else
    if size(txt,1)==1 && sizes(1)>1
        txt=txt';
    end
    txt(end+1:sizes(1),:)={''};
    txt(sizes(1)+1:end,:)=[];
    add=size(txt,2);
    txt=reshape(txt,sizes(1),add);
end
if nargin<4 || isempty(header)
    header={};
else
    sizh=size(header);
    if sizh(2)<sizes(2)+add
        toadd=min(add,sizes(2)+add-sizh(2));
        header=[repmat({''},[sizh(1),toadd]),header];
    end
    header(:,end+1:sizes(2)+add)={''};
    header(:,sizes(2)+add+1:end)=[];
    try
        header=reshape(header,sizh(1),sizes(2)+add);
    catch %old bug. Should not be reached any more
        header=reshape([{''},header],sizh(1),sizes(2)+add);
    end
end
if isempty(data)
    celldata=cell.empty(0,size(header,2)-size(txt,2));
elseif iscell(data)
    celldata=data;
elseif isnumeric(data)
    celldata=textscan([num2str(data),repmat(' ',sizes(1),1)]','%s');
    celldata=reshape(celldata{1},fliplr(sizes))'; % more portable, may be slightly slower
else
    error('luin:scrivi','scrivi may be used only with numeric or cell matrix');
end


% tic;
%     celldata1=reshape(strread([num2str(data),repmat('
%     ',sizes(1),1)]','%s'),fliplr(sizes))'; % old version
% toc;
% % tic;
% celldata= arrayfun(@num2str, data, 'UniformOutput', 0); %much slower
% toc;


if nargout>1
    celldata=[header;[txt, celldata]];
    byt=0;
else
    byt=writefile(fname,[header;[txt, celldata]],sep);
end
end