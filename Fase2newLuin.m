function Fase2newLuin(pn,doRadB,checkInCh,orientation)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
if nargin<4 || isempty(orientation)
    orientation='l'; %'left', 'right', 'up', 'down'
end
% BorderSize=10;
BorderSize=7; %bordo aggiunto per padding, usato per il fit gaussiano degli spot; pixels
MinimumLenght=0; %non usato, in implementazione precedente minima distanza percorsa in una traiettoria in micron
PixelSize=0.228; %um
min_n=15; %lunghezza minima traiettoria, frames
av_win=5; %dimensione media finestra mobile
min_asp_r=0; %somma spostamenti lungo x / somma spostamenti lungo y
dt=0.1; %sec
vmin=0.5; %vmin per considerare sottotraiettoria in moto (um/sec!)
vmint=0; %<|v|>min per accettare traiettoria (um/sec); se le vogliamo non ferme.
tc=4; %tempo minimo stop o go, frames
mindist=0; %pixels distanza minima dal bordo
colorcode={'red','black','green'};
maxstop=1.01; %massima frazione della traiettoria in cui la particella sta ferma.
appx='_tc2_minn15';
% BorderSize=7; %bordo aggiunto per padding, non usato ora; pixels
% MinimumLenght=10; %non usato, in implementazione precedente minima distanza percorsa in una traiettoria in micron
% PixelSize=0.228; %um
% min_n=15; %lunghezza minima traiettoria, frames
% av_win=5; %dimensione media finestra mobile
% min_asp_r=2.5; %somma spostamenti lungo x / somma spostamenti lungo y
% dt=0.106; %sec
% vmin=1.25; %(um/sec!)
% tc=4; %tempo minimo stop o go, frames
% mindist=3; %pixels distanza minima dal bordo
% channelwidth=34; %min, pixels
% colorcode={'red','black','green'};
% maxstop=0.8; %oscuro; c'entra col tempo totale di go e/o stop rispetto ai tempi di switching.
% appx='_l80_aw5tc4_v1_25';
if nargin<1 || isempty(pn)
%    olddir=cd('G:\Luos\2014\Research\Teresa\Esempio');
    olddir=pwd;
else
    olddir=cd(pn);
end
FileList=dir;
if nargin<2 || isempty(doRadB)
    doRadB=false;
end
if nargin<3 || isempty(checkInCh)
    checkInCh=false;
end
if ~checkInCh
    channelwidth=Inf;
else
    channelwidth=checkInCh;
end

tic;
for Ifile=3:size(FileList,1)
    FileName=FileList(Ifile).name;
    if length(FileName)>3 && strcmp(FileName(end-2:end),'mat') && strcmp(FileName(1:2),'Se') && FileList(Ifile).bytes>300 %&& ~exist(['exL_',FileName],'file')
        PartInfoC2=[];
        PartInfoC1=[];
        PartInfoC={PartInfoC1,PartInfoC2};
        extrC1=struct('t',[],'x',[],'y',[],'xmax',[],'xmin',[],'ymax',[],'ymin',[],'vx',[],'vy',[],'t_av',[],'orientation','','x_av',[],'y_av',[],'vx_av',[],'vy_av',[],'v_av',[],'SelId',[],'SelIdInd',[]);
        extrC2=extrC1;
        extrC={extrC1,extrC2};
        if doRadB
            Channel=load(['CH_',FileName]);
            Channel=Channel.Channel;
        end
        Tracking=load(FileName);
        Tracking=Tracking.Tracking;
        %% C1
        for chn=1:2
            if ~isempty(Tracking.(sprintf('objsC%d',chn)))
                Tracking.(sprintf('PartIdC%d',chn))=Tracking.(sprintf('objs_linkC%d',chn))(6,:);
                Tracking.(sprintf('SelIdC%d',chn))=[];
    %            Tracking.max=[];
    %            Tracking.min=[];
                II1=0;
                close all;
                ymax=max(Tracking.(sprintf('objs_linkC%d',chn))(2,:));
                ymin=min(Tracking.(sprintf('objs_linkC%d',chn))(2,:));
                delta=ymax-ymin-channelwidth;
                if checkInCh && delta>1 % controllo di stare dentro al canale, ipotizzando di avere i bordi del canale più intensi...
                    ymin=floor(ymin);
                    ymax=ceil(ymax);
                    delta=ceil(delta);
                    hnpts=hist(Tracking.(sprintf('objs_linkC%d',chn))(2,:),ymin:ymax);
                    [~, deltay]=max(hnpts(1:delta));
                    ymin=ymin+deltay-1;
                    [~, deltay]=max(hnpts(end:-1:end-delta));
                    ymax=ymax-deltay+1;
    %                 npnts=arrayfun(@(vec) sum((Tracking.(sprintf('objs_linkC%d',chn))(2,:)>vec)&(Tracking.(sprintf('objs_linkC%d',chn))(2,:)<vec+channelwidth+1)),floor(ymin):(ceil(ymax)-channelwidth));
    %                 [~,ymin]=max(npnts);
    %                 ymax=ymin+channelwidth+1;
                end
                
                for II=unique(Tracking.(sprintf('PartIdC%d',chn)))
                    ind=find(Tracking.(sprintf('PartIdC%d',chn))==II);
                    nframes=max(Tracking.(sprintf('objs_linkC%d',chn))(5,:));
                    if length(ind)>min_n
                        ymean=mean(Tracking.(sprintf('objs_linkC%d',chn))(2,ind));
                        if ymean>ymin+mindist && ymean<ymax-mindist
                            x_av=movingmean(Tracking.(sprintf('objs_linkC%d',chn))(1,ind),av_win);
                            y_av=movingmean(Tracking.(sprintf('objs_linkC%d',chn))(2,ind),av_win);
                            t_av=movingmean(Tracking.(sprintf('objs_linkC%d',chn))(5,ind),av_win);
                            
                            vx_av=diff(x_av)./diff(t_av);
                            vx_avav=mean(abs(vx_av));
                            vy_av=diff(y_av)./diff(t_av);
                            vy_avav=mean(abs(vy_av));
                            v_av=sqrt(vx_av.^2+vy_av.^2);
                            v_avav=mean(v_av);
                            switch(orientation(1))
                                case 'l'
                                    orn=-sign(vx_av);
                                    cond_angle=vx_avav>min_asp_r*vy_avav;
                                case 'u'
                                    orn=-sign(vy_av);
                                    cond_angle=vy_avav>min_asp_r*vx_avav;
                                case 'd'
                                    cond_angle=vy_avav>min_asp_r*vx_avav;
                                    orn=sign(vy_av);
                                otherwise
                                    orn=sign(vx_av);
                                    cond_angle=vx_avav>min_asp_r*vy_avav;
                            end
                            v_av=v_av.*orn;
                            
                            [S,starts,nsub]=StopFwRv(v_av,vmin*dt/PixelSize,tc,t_av);
                            de_t=diff(Tracking.(sprintf('objs_linkC%d',chn))(5,ind));
                            if v_avav>vmint*dt/PixelSize && cond_angle && sum(de_t(~S))<maxstop*sum(de_t) %parametro diverso per vmin?
                                II1=II1+1;
                                Tracking.(sprintf('SelIdC%d',chn))(II1)=II;
                                extrC{chn}(II1).SelId=II;
                                extrC{chn}(II1).SelIdInd=ind;
                                extrC{chn}(II1).x=Tracking.(sprintf('objs_linkC%d',chn))(1,ind);
                                extrC{chn}(II1).y=Tracking.(sprintf('objs_linkC%d',chn))(2,ind);
                                extrC{chn}(II1).t=Tracking.(sprintf('objs_linkC%d',chn))(5,ind);
                                
                                [extrC{chn}(II1).xmin(1), extrC{chn}(II1).xmin(2)]=min(extrC{chn}(II1).x);
                                [extrC{chn}(II1).xmax(1), extrC{chn}(II1).xmax(2)]=max(extrC{chn}(II1).x);
                                [extrC{chn}(II1).ymin(1), extrC{chn}(II1).ymin(2)]=min(extrC{chn}(II1).y);
                                [extrC{chn}(II1).ymax(1), extrC{chn}(II1).ymax(2)]=max(extrC{chn}(II1).y);
    %                            extrC{chn}(II1).min(2)=extrC{chn}(II1).t(b);
    %                            extrC{chn}(II1).max(2)=ind(b);
                                extrC{chn}(II1).vx=diff(extrC{chn}(II1).x)./diff(extrC{chn}(II1).t);
                                extrC{chn}(II1).vy=diff(extrC{chn}(II1).y)./diff(extrC{chn}(II1).t);
                                extrC{chn}(II1).orientation=orientation;
                                extrC{chn}(II1).t_av=t_av;
                                extrC{chn}(II1).x_av=x_av;
                                extrC{chn}(II1).y_av=y_av;
                                extrC{chn}(II1).vx_av=vx_av;
                                extrC{chn}(II1).vy_av=vy_av;
                                extrC{chn}(II1).v_av=v_av;
                                extrC{chn}(II1).S=S;
                                extrC{chn}(II1).starts=starts;
                                extrC{chn}(II1).nsub=nsub;
                                plot(Tracking.(sprintf('objs_linkC%d',chn))(1,ind),Tracking.(sprintf('objs_linkC%d',chn))(2,ind));axis image;drawnow
                                hold on;
                                for III=1:nsub
                                    plot(x_av(starts(III):starts(III+1)),y_av(starts(III):starts(III+1)),colorcode{S(starts(III))+2});
                                end
                                drawnow;
                                %                            plot(x_av,y_av,'red');drawnow;
                                
                            end
                            %{
    %                     vx=diff()./diff(Tracking.(sprintf('objs_linkC%d',chn))(5,ind));
    %                     vxav=mean(abs(vx));
    %                     if vxav>MinimumLenght/(PixelSize*nframes)
    %                         II1=II1+1;
    %                         Tracking.(sprintf('SelIdC%d',chn))(II1)=II;
    %                         [Tracking.min(1,II1) b]=min(Tracking.(sprintf('objs_linkC%d',chn))(1,ind));
    %                         Tracking.min(2,II1)=ind(b);
    %                         [Tracking.max(1,II1) b]=max(Tracking.(sprintf('objs_linkC%d',chn))(1,ind));
    %                         Tracking.max(2,II1)=ind(b);
    %                         plot(Tracking.(sprintf('objs_linkC%d',chn))(1,Tracking.(sprintf('PartIdC%d',chn))==II),Tracking.(sprintf('objs_linkC%d',chn))(2,Tracking.(sprintf('PartIdC%d',chn))==II));axis image;drawnow
    %                         hold on
    %                     end
                            %}
                        end
                    end
                end
                
                if numel(Tracking.(sprintf('SelIdC%d',chn)))>0 % isempty
                    if doRadB
                        Channel.(sprintf('C%d',chn))=padarray(Channel.(sprintf('C%d',chn)),[BorderSize BorderSize]);
                    end
                    print(1,['CH' int2str(chn) 'LL_',FileName(1:end-4),appx,'.png'],'-dpng')
                    %                close all
                    PartInfoC{chn}=cell(numel(Tracking.(sprintf('SelIdC%d',chn))),5);
                    I1=0;
                    if doRadB
                        for I=Tracking.(sprintf('SelIdC%d',chn))
                            I1=1+I1;
                            PartX=Tracking.(sprintf('objs_linkC%d',chn))(1,Tracking.(sprintf('PartIdC%d',chn))==I);
                            PartY=Tracking.(sprintf('objs_linkC%d',chn))(2,Tracking.(sprintf('PartIdC%d',chn))==I);
                            T=Tracking.(sprintf('objs_linkC%d',chn))(5,Tracking.(sprintf('PartIdC%d',chn))==I);
                            Part=zeros(BorderSize*2+1,BorderSize*2+1,numel(T));
                            II1=0;
                            for II=T
                                II1=II1+1;
                                Part(:,:,II1)=Channel.(sprintf('C%d',chn))(floor(PartY(II1))+1:floor(PartY(II1))+1+2*BorderSize,floor(PartX(II1))+1:floor(PartX(II1))+1+2*BorderSize,II);
                                
                            end
                            if size(Part,3)>20
                                [a] = GaussFit_fixed(MovAv(Part,10),1,1+floor(size(Part,1)/2),1+floor(size(Part,1)/2),[max(Part(:))-min(Part(:)) min(Part(:)) 2]);
                            else
                                [a] = GaussFit_fixed(mean(Part,3),1,1+floor(size(Part,1)/2),1+floor(size(Part,1)/2),[max(Part(:))-min(Part(:)) min(Part(:)) 2]);
                            end
                            %                 [a] = GaussFit_Teresa(Part,PartX-floor(PartX)+1+BorderSize,PartY-floor(PartY)+1+BorderSize);
                            PartInfoC{chn}{I1,1}=PartX;
                            PartInfoC{chn}{I1,2}=PartY;
                            PartInfoC{chn}{I1,3}=T;
                            PartInfoC{chn}{I1,4}=Part;
                            PartInfoC{chn}{I1,5}=a;
                        end
                    end
                end
            end
        end
        extrC1=extrC{1}; %#ok<NASGU>
        extrC2=extrC{2}; %#ok<NASGU>
        PartInfoC1=PartInfoC{1}; %#ok<NASGU>
        PartInfoC2=PartInfoC{2}; %#ok<NASGU>
        save(['exL_',FileName(1:end-4),appx,'.mat'],'PartInfoC2','PartInfoC1','extrC1','extrC2')
    end
end
toc;
cd(olddir);
%%
% %% C2
%         if ~isempty(Tracking.objsC2)
%             Tracking.PartIdC2=Tracking.objs_linkC2(6,:);
%             Tracking.SelIdC2=[];
% %            Tracking.max=[];
% %            Tracking.min=[];
%             II1=0;
%             close all;
%             ymax=max(Tracking.objs_linkC2(2,:));
%             ymin=min(Tracking.objs_linkC2(2,:));
%             delta=ymax-ymin-channelwidth;
%             if delta>1
%                 ymin=floor(ymin);
%                 ymax=ceil(ymax);
%                 delta=ceil(delta);
%                 hnpts=hist(Tracking.objs_linkC2(2,:),ymin:ymax);
%                 [~, deltay]=max(hnpts(1:delta));
%                 ymin=ymin+deltay-1;
%                 [~, deltay]=max(hnpts(end:-1:end-delta));
%                 ymax=ymax-deltay+1;
% %                 npnts=arrayfun(@(vec) sum((Tracking.objs_linkC2(2,:)>vec)&(Tracking.objs_linkC2(2,:)<vec+channelwidth+1)),floor(ymin):(ceil(ymax)-channelwidth));
% %                 [~,ymin]=max(npnts);
% %                 ymax=ymin+channelwidth+1;
%             end
%
%             for II=unique(Tracking.PartIdC2)
%                 ind=find(Tracking.PartIdC2==II);
%                 nframes=max(Tracking.objs_linkC2(5,:));
%                 if length(ind)>min_n
%                     ymean=mean(Tracking.objs_linkC2(2,ind));
%                     if ymean>ymin+mindist && ymean<ymax-mindist
%                         x_av=movingmean(Tracking.objs_linkC2(1,ind),av_win);
%                         y_av=movingmean(Tracking.objs_linkC2(2,ind),av_win);
%                         t_av=movingmean(Tracking.objs_linkC2(5,ind),av_win);
%
%                         vx_av=diff(x_av)./diff(t_av);
%                         vx_avav=mean(abs(vx_av));
%                         vy_av=diff(y_av)./diff(t_av);
%                         vy_avav=mean(abs(vy_av));
%
%                         [S,starts,nsub]=StopFwRv(vx_av,vmin*dt/PixelSize,tc,t_av);
%                         de_t=diff(Tracking.objs_linkC2(5,ind));
%                         if vx_avav>vmint*dt/PixelSize && vx_avav>min_asp_r*vy_avav && sum(de_t(~S))<maxstop*sum(de_t)
%                             II1=II1+1;
%                             Tracking.SelIdC2(II1)=II;
%                             extrC2(II1).SelId(II1)=II;
%                             extrC2(II1).SelIdInd=ind;
%                             extrC2(II1).x=Tracking.objs_linkC2(1,ind);
%                             extrC2(II1).y=Tracking.objs_linkC2(2,ind);
%                             extrC2(II1).t=Tracking.objs_linkC2(5,ind);
%
%                             [extrC2(II1).xmin(1) extrC2(II1).xmin(2)]=min(extrC2(II1).x);
% %                            extrC2(II1).min(2)=extrC2(II1).t(b);
%                             [extrC2(II1).xmax(1) extrC2(II1).xmax(2)]=max(extrC2(II1).x);
% %                            extrC2(II1).max(2)=ind(b);
%                             extrC2(II1).vx=diff(extrC2(II1).x)./diff(extrC2(II1).t);
%                             extrC2(II1).vy=diff(extrC2(II1).y)./diff(extrC2(II1).t);
%                             extrC2(II1).t_av=t_av;
%                             extrC2(II1).x_av=x_av;
%                             extrC2(II1).y_av=y_av;
%                             extrC2(II1).vx_av=vx_av;
%                             extrC2(II1).vy_av=vy_av;
%                             extrC2(II1).S=S;
%                             extrC2(II1).starts=starts;
%                             extrC2(II1).nsub=nsub;
%                             plot(Tracking.objs_linkC2(1,ind),Tracking.objs_linkC2(2,ind));axis image;drawnow
%                             hold on;
%                             for III=1:nsub
%                                 plot(x_av(starts(III):starts(III+1)),y_av(starts(III):starts(III+1)),colorcode{S(starts(III))+2});
%                             end
% %                            plot(x_av,y_av,'red');drawnow;
%
%                         end
%                         %{
%     %                     vx=diff()./diff(Tracking.objs_linkC2(5,ind));
%     %                     vxav=mean(abs(vx));
%     %                     if vxav>MinimumLenght/(PixelSize*nframes)
%     %                         II1=II1+1;
%     %                         Tracking.SelIdC2(II1)=II;
%     %                         [Tracking.min(1,II1) b]=min(Tracking.objs_linkC2(1,ind));
%     %                         Tracking.min(2,II1)=ind(b);
%     %                         [Tracking.max(1,II1) b]=max(Tracking.objs_linkC2(1,ind));
%     %                         Tracking.max(2,II1)=ind(b);
%     %                         plot(Tracking.objs_linkC2(1,Tracking.PartIdC2==II),Tracking.objs_linkC2(2,Tracking.PartIdC2==II));axis image;drawnow
%     %                         hold on
%     %                     end
%                         %}
%                     end
%                 end
%             end
%
%             if numel(Tracking.SelIdC2)>0 % >1??? isempty?
%                 if doRadB
%                     Channel.C2=padarray(Channel.C2,[BorderSize BorderSize]);
%                 end
%                 print(1,['CH2LL_',FileName(1:end-4),appx,'.png'],'-dpng')
%                 %                close all
%                 PartInfoC2=cell(numel(Tracking.SelIdC2),5);
%                 I1=0;
%                 if doRadB
%
%                     for I=Tracking.SelIdC2
%                         I1=1+I1;
%                         PartX=Tracking.objs_linkC2(1,Tracking.PartIdC2==I);
%                         PartY=Tracking.objs_linkC2(2,Tracking.PartIdC2==I);
%                         T=Tracking.objs_linkC2(5,Tracking.PartIdC2==I);
%                         Part=zeros(BorderSize*2+1,BorderSize*2+1,numel(T));
%                         II1=0;
%                         for II=T
%                             II1=II1+1;
%                             Part(:,:,II1)=Channel.C2(floor(PartY(II1))+1:floor(PartY(II1))+1+2*BorderSize,floor(PartX(II1))+1:floor(PartX(II1))+1+2*BorderSize,II);
%
%                         end
%                         if size(Part,3)>20
%                             [a] = GaussFit_fixed(MovAv(Part,10),1,1+floor(size(Part,1)/2),1+floor(size(Part,1)/2),[max(Part(:))-min(Part(:)) min(Part(:)) 2]);
%                         else
%                             [a] = GaussFit_fixed(mean(Part,3),1,1+floor(size(Part,1)/2),1+floor(size(Part,1)/2),[max(Part(:))-min(Part(:)) min(Part(:)) 2]);
%                         end
%                         %                 [a] = GaussFit_Teresa(Part,PartX-floor(PartX)+1+BorderSize,PartY-floor(PartY)+1+BorderSize);
%                         PartInfoC2{I1,1}=PartX;
%                         PartInfoC2{I1,2}=PartY;
%                         PartInfoC2{I1,3}=T;
%                         PartInfoC2{I1,4}=Part;
%                         PartInfoC2{I1,5}=a;
%                     end
%                 end
%             end
%         end
%
%{
%         if ~isempty(Tracking.objsC2)
%         Tracking.PartIdC2=Tracking.objs_linkC2(6,:);
%         Tracking.SelIdC2=[];
%         II1=0;
%         close all
%         for II=unique(Tracking.PartIdC2)
%
%             if max(Tracking.objs_linkC2(1,Tracking.PartIdC2==II))-min(Tracking.objs_linkC2(1,Tracking.PartIdC2==II))>MinimumLenght/PixelSize
%                 II1=II1+1;
%                 Tracking.SelIdC2(II1)=II;
% %                 1
%                 plot(Tracking.objs_linkC2(1,Tracking.PartIdC2==II),Tracking.objs_linkC2(2,Tracking.PartIdC2==II));axis image;drawnow
%                 hold on
%             end
%
%         end
%         if numel(Tracking.SelIdC2)>1
%             Channel.C2=padarray(Channel.C2,[BorderSize BorderSize]);
%
%             print(1,['CH2_',FileName(1:end-3),'.png'],'-dpng')
%             close all
%             PartInfoC2=cell(numel(Tracking.SelIdC2),5);
%             I1=0;
%             for I=Tracking.SelIdC2
%                 I1=1+I1;
%                 PartX=Tracking.objs_linkC2(1,Tracking.PartIdC2==I);
%                 PartY=Tracking.objs_linkC2(2,Tracking.PartIdC2==I);
%                 T=Tracking.objs_linkC2(5,Tracking.PartIdC2==I);
%                 Part=zeros(BorderSize*2+1,BorderSize*2+1,numel(T));
%                 II1=0;
%                 for II=T
%                     II1=II1+1;
%                     Part(:,:,II1)=Channel.C2(floor(PartY(II1))+1:floor(PartY(II1))+1+2*BorderSize,floor(PartX(II1))+1:floor(PartX(II1))+1+2*BorderSize,II);
%                 end
%                 if size(Part,3)>20
%                     [a] = GaussFit_fixed(MovAv(Part,10),1,1+floor(size(Part,1)/2),1+floor(size(Part,1)/2),[max(Part(:))-min(Part(:)) min(Part(:)) 2]);
%                 else
%                     [a] = GaussFit_fixed(mean(Part,3),1,1+floor(size(Part,1)/2),1+floor(size(Part,1)/2),[max(Part(:))-min(Part(:)) min(Part(:)) 2]);
%                 end
%                 %                 [a] = GaussFit_Teresa(Part,PartX-floor(PartX)+1+BorderSize,PartY-floor(PartY)+1+BorderSize);
%                 PartInfoC2{I1,1}=PartX;
%                 PartInfoC2{I1,2}=PartY;
%                 PartInfoC2{I1,3}=T;
%                 PartInfoC2{I1,4}=Part;
%                 PartInfoC2{I1,5}=a;
%             end
%         end
%         end
%}
