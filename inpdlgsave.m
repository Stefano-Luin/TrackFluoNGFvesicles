function [out,answ,varstr]=inpdlgsave(namef,varstr,defs,typestr,texts,title,prefb,optdlg,lines)
% out=inpdlgsave(namef,varstr,defs,typestr,texts,title,prefb,optdlg,lines);
% prefb=XY_2, X=out in pref, Y=answ in pref;
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)

if nargin<9||isempty(lines)
    lines=1;
end
if nargin<2
    error('Luin:inpdlgsave','inpdlgsave called with too few arguments');
end

if isempty(varstr)
    if nargin>4 && ~isempty(texts)
        if ischar(texts)
            texts=cellstr(texts);
        end
        varstr = texts;
    else
        error('Luin:inpdlgsave','inpdlgsave called with too few arguments');
    end
end

if ischar(varstr)
    varstr=cellstr(varstr);
end

if nargin<8 || isempty(optdlg)
    optdlg.Resize='on';
    optdlg.WindowStyle='normal';
    optdlg.Interpreter='none';
end

[pdir,pf,ext]=fileparts(namef);
if ~strcmp(ext,'.mat')
    namef=fullfile(pdir,[pf '.mat']);
end

sv=size(varstr);
if nargin<3 ||isempty(defs)
    defs=repmat({''},sv);
end

if nargin<4||isempty(typestr)
    typestr=repmat('s',sv);
end
if nargin<5 || isempty(texts)
    texts=varstr;
end
if nargin<7||isempty(prefb)
    prefb=false;
end
if nargin<6||isempty(title)
    title=pf;
end
answ={};
if mod(prefb,2) && pref('is','inpdlg',pf)
    answ=pref('get','inpdlg',pf);
% elseif exist(fullfile(pwd,namef),'file')
%     answ=load(fullfile(pwd,namef));
elseif ~isempty(ls(namef))
    answ=load(namef);
    if isfield(answ,'answ')
        answ=answ.answ;
    else
        answ={};
    end
end
if size(answ)==sv
    defs=answ;
end
cl=find(typestr=='c');
for in=cl
    if iscell(defs{in})
        defs{in}(cellfun(@isnumeric,defs{in}))=cellfun(@num2str,defs{in}(cellfun(@isnumeric,defs{in})),'UniformOutput',false);
        if ~isempty(defs{in})
            temp='{';
            for i=1:size(defs{in},1)
                temp=[temp,sprintf('''%s'',',defs{in}{i,:}),';']; %#ok<AGROW>
                temp(end-1)=[];
            end
            temp(end)='}';
            defs{in}=temp;
        else
            defs{in}='';
        end
    end
end
answ=inputdlg(texts, title, lines, defs,optdlg);
answ=reshape(answ,sv);
drawnow; pause(0.01);

if verLessThan('matlab','8.4')
    varstr=genvarname(varstr); 
else
    varstr=matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(varstr));
end

for c=1:prod(sv)
    switch typestr(c)
        case 'n'
            out.(varstr{c})=str2num(answ{c}); %#ok<ST2NM>
        case 'd'
            out.(varstr{c})=str2double(answ{c});
        case 'c' %cell
            if isempty(answ{c})
                out.(varstr{c})={};
            else
                out.(varstr{c})=eval(answ{c});
            end
            answ{c}=out.(varstr{c});
        case {'l','b'}
            out.(varstr{c})=logical(str2num(answ{c})); %#ok<ST2NM>
        otherwise
            out.(varstr{c})=answ{c};
    end
end

if mod(prefb,2)
    pref('set', 'inpdlg',pf,answ);
end
try %#ok<TRYNC>
    save(checkdir(namef),'answ');
end
if prefb>=2
    pref('set','Luos',pf,out);
end
