% load('info_l80_aw5tc4_v1_25.mat');
% Copyright (C) 2012 - 2022 Domenica Convertino, Stefano Luin (s.luin@sns.it)
llls= 'exp';
par='_aw5tc4_v05';

makehist([[extr.v_avsn,-[extr.v_avdx]]',[sqrt([extr.v_avsn_s2s]),sqrt([extr.v_avdx_s2s])]',[extr.dtsn,extr.dtdx]'],(-10:0.25:10),[llls par],'v_go_s',1);
makehist([[extr.v_avsn,-[extr.v_avdx]]',[sqrt([extr.v_avsn_s2m]),sqrt([extr.v_avdx_s2m])]',[extr.dtsn,extr.dtdx]'],(-10:0.25:10),[llls par],'v_go_m',1);
%makehist([[extr.v_avdx]',sqrt([extr.v_avdx_s2s])',[extr.dndx]'],(-10:0.25:10),[llls par],'vdx',1);
%makehist([[extr.v_avsn]',sqrt([extr.v_avsn_s2s])',[extr.dnsn]'],(-10:0.25:10),[llls par],'vsn',1);
%makehist([[extr.v_aval]',sqrt([extr.v_aval_s2s])',[extr.dtal]'],(-10:0.25:10),[llls par],'val_s',1);
%makehist([[extr.v_aval]',sqrt([extr.v_aval_s2m])',[extr.dtal]'],(-10:0.25:10),[llls par],'val_m',1);
%makehist([-[extr.vmed]',[extr.vmed_s]',[extr.t_tot]'],(-10:0.25:10),[llls par],'vmed',1);
%makehist([-[extr.vmax]',zeros(size([extr.vmax]')),[extr.t_tot]'],(-10:0.25:10),[llls par],'vmax',1);
%makehist([-[extr.v_old]',zeros(size([extr.v_old]')),[extr.n_tot]'],(-10:0.25:10),[llls par],'vtot',1);

