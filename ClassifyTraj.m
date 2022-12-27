% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
OutName='outb';
vminosc=0.1; %<v>min per traiettoria non oscillante.
vmint=0.03; %vmin per immobile, di controllo.
maxstop=0.9; %massima frazione della traiettoria in cui la particella sta ferma.
ndiffmax=2;
rationdiff=3;
minlength=30; %frame
minshiftmob=5; %pixel
minshiftosc=10; %pixel
extr(end).class='';
for ii=1:length(extr)
    if extr(ii).n_tot<minlength
        extr(ii).class='short';
    elseif ((1-maxstop)*extr(ii).dtal>=extr(ii).dtsn+extr(ii).dtdx && ...
            abs(extr(ii).vmed)<vmint) || ... %nota: dt è solitamente maggiore della somma dei dt, essendo questi ultimi calcolati su t_av.
            ((extr(ii).x_av(end)-extr(ii).x_av(1))^2+(extr(ii).y_av(end)-extr(ii).y_av(1))^2<minshiftmob^2 || ...
            abs(extr(ii).ornv*[extr(ii).x_av(end)-extr(ii).x_av(1);extr(ii).y_av(end)-extr(ii).y_av(1)])<minshiftmob/sqrt(2))
        extr(ii).class='immobile';
    elseif (extr(ii).dnsn>0 && extr(ii).dndx>0 && ...
            abs(sum(extr(ii).S(extr(ii).starts(1:end-1))))<=max(ndiffmax,length(extr(ii).starts)/rationdiff)...
            && abs(extr(ii).vmed)<vminosc)||...
            ((extr(ii).x_av(end)-extr(ii).x_av(1))^2+(extr(ii).y_av(end)-extr(ii).y_av(1))^2<minshiftosc^2 || ...
            abs(extr(ii).ornv*[extr(ii).x_av(end)-extr(ii).x_av(1);extr(ii).y_av(end)-extr(ii).y_av(1)])<minshiftosc/sqrt(2))
        extr(ii).class='oscillating';
    elseif extr(ii).v_old>0
        extr(ii).class='dx';
    else
        extr(ii).class='sn';
    end
end

shr=arrayfun(@(s) strcmp(s.class,'short'),extr);
osc=arrayfun(@(s) strcmp(s.class,'oscillating'),extr);
dx=arrayfun(@(s) strcmp(s.class,'dx'),extr);
sn=arrayfun(@(s) strcmp(s.class,'sn'),extr);
imbl=arrayfun(@(s) strcmp(s.class,'immobile'),extr);
%%
classes={'short','osc','dx','sn','imbl'};
classId={shr,osc,dx,sn,imbl};
t=struct('n_tot',[],'t_tot',[],'shift',[],'n_traj',[]);
t(length(classes))=t;
for jj=1:length(classes)
    pplot=find(classId{jj});
    Fase4(extr(pplot),OutName,classes{jj});
    
    f1= figure('Name',[classes{jj} '_xy']); hold on;
    axis image;
    for ii=pplot
        plot(extr(ii).x,extr(ii).y);
    end
    print(f1,[OutName classes{jj} '_xy'],'-dpng');
    
    f2=figure('Name',[classes{jj} '_tx']); hold on;
    for ii=pplot
        extr(ii).ds=[extr(ii).vxmax./extr(ii).vmax,extr(ii).vymax./extr(ii).vmax]*[extr(ii).x-extr(ii).x(1);extr(ii).y-extr(ii).y(1)];
        plot(extr(ii).t-extr(ii).t(1),extr(ii).ds);
%         plot(extr(ii).t-extr(ii).t(1),extr(ii).x-extr(ii).x(1));
    end
    print(f2,[OutName classes{jj} '_tx'],'-dpng');
    
    t(jj).n_tot=sum([extr(pplot).n_tot]);
    t(jj).t_tot=sum([extr(pplot).t_tot]);
    t(jj).shift=[extr(pplot).t_tot]*[extr(pplot).v_old]';
    t(jj).n_traj=length(pplot);
end
tc=struct2cell(t);
tc=cellfun(@replaceempty,tc,'UniformOutput',false);
tc=cell2mat(squeeze(tc)');
scrivi([OutName 'ClassTraj.dat'],tc,classes',fieldnames(t)');
save([OutName 'extr.mat'],'extr','classes','classId','t');

function oo=replaceempty(ii)
if isempty(ii)
    oo=NaN;
else
    oo=ii;
end
end