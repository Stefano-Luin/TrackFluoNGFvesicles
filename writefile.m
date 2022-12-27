function output=writefile(filename,cells,sep)
%output=writefile(filename,cells)
%writefile save the first 2D cell array
%   see also APRI SCRIVI READFILE OPENASCII
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
if nargin<3||isempty(sep)
    sep='\t';
%    sep=', ';
%    filename=[filename '.csv'];
end
fid=openw(filename,'wt');
sizes=size(cells);
if length(sizes)~=2
    warning('Luin:writefile','writefile save the first 2D cell array');
end
formatstring=repmat(['%s' sep],1,sizes(2)-1);
formatstring=[formatstring '%s\n'];
outp=zeros(1,sizes(1));
for row=1:sizes(1)
    outp(row)=fprintf(fid,formatstring,cells{row,:,1});
end
if nargout>0, output=sum(outp); end
fclose(fid);
end