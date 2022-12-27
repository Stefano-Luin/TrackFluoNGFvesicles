function [nameout,done]=checkdir(nameout)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
ddir=fileparts(nameout);
done=false;
if ~isempty(ddir)&&~exist(ddir,'dir')
    mkdir(ddir);
    done=true;
end
