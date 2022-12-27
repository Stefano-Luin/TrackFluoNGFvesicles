function out=ffpath(in)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
out=fileparts(in);
global usrpth
if isempty(usrpth)
    if pref('is','Luos','usrpth')
        usrpth=pref('get','Luos','usrpth');
    else
        usrpth=userpath;
        usrpth(end)='';
    end
end
switch in
    case 'loci_tools.jar'
        try %#ok<TRYNC>
            out=fullfile(usrpth,'Carmine\Mfile\Suppl'); %to be changed
        end
end
% out='D:\Download';
