function parDom=AskDomenica()
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
% [tmpobj] = fo5_rp(double(Channel.C1(:,:,j))-SmoothGauss(Channel.C1(:,:,j),100), 'spatialfilter', [10 10], 0.99, 'radial');
% Tracking.objs_linkC1 = nnlink_rp(Tracking.objsC1, 30, 2);
clear varstr defs texts typestr
c=0;
c=c+1;
varstr{c}='GSmooth';
defs{c}='100';
texts{c}='Gauss smoothing for background subtract';
typestr(c)='d';
c=c+1;
varstr{c}='processopt';
defs{c}='spatialfilter';
texts{c}='option for processing to find neighborhoods';
typestr(c)='s';
c=c+1;
varstr{c}='processparam';
defs{c}='[10,10]';
texts{c}='Parameters for processing';
typestr(c)='n';
c=c+1;
varstr{c}='thresh';
defs{c}='0.99';
texts{c}='intensity threshold';
typestr(c)='d';
c=c+1;
varstr{c}='fitstr';
defs{c}='radial';
texts{c}='String (or array of strings) that selects the fitting option';
typestr(c)='s';
c=c+1;
varstr{c}='step';
defs{c}='30';
texts{c}='step: cull all links that are greater than sqrt(step*memory) away';
typestr(c)='d';
c=c+1;
varstr{c}='memory';
defs{c}='2';
texts{c}='memory: attempt to link particles losing them';
typestr(c)='d';

parDom=inpdlgsave('parDom',varstr,defs,typestr,texts,'Localization and linking paramters',3);