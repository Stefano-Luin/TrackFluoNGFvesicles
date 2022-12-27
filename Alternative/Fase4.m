function [r,histtype] = Fase4(extr,llls,par,h_todo)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
histtype={'v_go_s','v_go_m','vdx','vsn','val_s','val_m','vmed','vmax','vtot'};
if nargin<4 || isempty(h_todo) || ~iscell(h_todo)
    h_todo=histtype;
end
if strcmp(h_todo{1},'not')
    h_todo=setdiff(histtype,h_todo);
end
a=ismember(h_todo,histtype);
if ~all(a)
    warning('luin:Fase4',[sprintf('%s ',h_todo{~a}) 'not valid histogram type(s)']);
end
if nargin<2 || isempty(llls)
    llls= 'exp';
end
if nargin<3 || ~ischar(par)
    par='_aw5tc4_v05';
end
if nargin<1 || isempty(extr)
    try
    extr=evalin('caller','extr');
    catch
        error('Luin:Fase4','Fase4 called with too few arguments');
    end
end
c=0;
if ismember('v_go_s',h_todo)
    c=c+1;
    r{c}=makehist([[extr.v_avsn,-[extr.v_avdx]]',[sqrt([extr.v_avsn_s2s]),sqrt([extr.v_avdx_s2s])]',[extr.dtsn,extr.dtdx]'],(-10:0.25:10),[llls par],'v_go_s',1);
end
if ismember('v_go_m',h_todo)
    c=c+1;
    r{c}=makehist([[extr.v_avsn,-[extr.v_avdx]]',[sqrt([extr.v_avsn_s2m]),sqrt([extr.v_avdx_s2m])]',[extr.dtsn,extr.dtdx]'],(-10:0.25:10),[llls par],'v_go_m',1);
end
if ismember('vdx',h_todo)
    c=c+1;
    r{c}=makehist([[extr.v_avdx]',sqrt([extr.v_avdx_s2s])',[extr.dndx]'],(-10:0.25:10),[llls par],'vdx',1);
end
if ismember('vsn',h_todo)
    c=c+1;
    r{c}=makehist([[extr.v_avsn]',sqrt([extr.v_avsn_s2s])',[extr.dnsn]'],(-10:0.25:10),[llls par],'vsn',1);
end
if ismember('val_s',h_todo)
    c=c+1;
    r{c}=makehist([[extr.v_aval]',sqrt([extr.v_aval_s2s])',[extr.dtal]'],(-10:0.25:10),[llls par],'val_s',1);
end
if ismember('val_m',h_todo)
    c=c+1;
    r{c}=makehist([[extr.v_aval]',sqrt([extr.v_aval_s2m])',[extr.dtal]'],(-10:0.25:10),[llls par],'val_m',1);
end
if ismember('vmed',h_todo)
    c=c+1;
    r{c}=makehist([-[extr.vmed]',[extr.vmed_s]',[extr.t_tot]'],(-10:0.25:10),[llls par],'vmed',1);
end
if ismember('vmax',h_todo)
    c=c+1;
    r{c}=makehist([-[extr.vmax]',zeros(size([extr.vmax]')),[extr.t_tot]'],(-10:0.25:10),[llls par],'vmax',1);
end
if ismember('vtot',h_todo)
    c=c+1;
    r{c}=makehist([-[extr.v_old]',zeros(size([extr.v_old]')),[extr.n_tot]'],(-10:0.25:10),[llls par],'vtot',1);
end

