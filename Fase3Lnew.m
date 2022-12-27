function Fase3Lnew(FolderList,appx,exL,calcexc,tip,InfoName)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
if nargin<5 || isempty(tip), tip=0; end
if nargin<4 || isempty(calcexc), calcexc=0; end
if nargin<6 || isempty(InfoName)
    InfoName='Info';
end
av_win=5;
vmin=0.5; %vmin per considerare sottotraiettoria in moto (um/sec)
tc=4;

%BorderSize=7;
if nargin<2 || isempty(appx) || ~ischar(appx)
    appx='';
end
if nargin<3 || isempty(exL) || ~ischar(exL)
    exL='exL_';
end
lappx=numel(appx);
lexL=numel(exL);
PixelSize=0.228;
A=importdata('Time.txt');
SeriesPos=A.data(:,1);
TimeFromAd=A.data(:,2);
% AllSeriesInd=SeriesInfo(:,1);
% SeriesAcqTime=SeriesInfo(:,2);
% SeriesTime=SeriesInfo(:,3);
% SeriesGain=SeriesInfo(:,4);
% SeriesEMGain=SeriesInfo(:,5);
% SeriesDelayTime=SeriesInfo(:,6);

RootFolder=pwd;
if tip, cd('Tip'); end
if nargin<1 || isempty(FolderList), FolderList=dir; end %da controllare queste righe se ~isempty(FolderList).

Itot=0;
%ItotR=0;
extr=[];
info={};
stat=warning('off','LUIN:wmean1');
warning('off','Luin:wmean');

for Ifolder=1:size(FolderList,1)
    if FolderList(Ifolder).isdir && exist(fullfile(RootFolder,[FolderList(Ifolder).name,'.lif']),'file')
        cd(FolderList(Ifolder).name)
        disp(['Fase3 in ' FolderList(Ifolder).name ' ' datestr(now)]);
        FileList=dir;
        for Ifile=3:size(FileList,1)
            FileName=FileList(Ifile).name;
            if length(FileName)>4+lappx+lexL && strcmp(FileName(1:lexL),exL) && strcmp(FileName((end-3-lappx):end),[appx '.mat'])
                disp(['Fase3: ' FileName ' ' datestr(now)]);
                %                SeriesNum=str2num(FileName(13:end-4));
                PartInfoC2=[];
                PartInfoC1=[];
                extrC1=[];
                extrC2=[];
                load(FileName); %da assegnare ad una variabile
                PartInfoC={PartInfoC1,PartInfoC2};
                extrC={extrC1,extrC2};
                %% C1 & C2
                for ii=1:2
                    if ~isempty(PartInfoC{ii})
                        if isempty(extrC{ii}), calc_exC=1;
                        else calc_exC=0;
                        end
                        SeriesNum=cell2mat(textscan(FileName, [exL 'Series_%d8' appx '.mat']));
                        try
                            SeriesInfo=ReadInfo(fullfile(RootFolder,[FolderList(Ifolder).name,'.lif']),SeriesNum,SeriesNum);
                        catch ME
                            warning('Luin:Fase3Lnew_java', 'Error in Java; retrying');
                            SeriesInfo=ReadInfo(fullfile(RootFolder,[FolderList(Ifolder).name,'.lif']),SeriesNum,SeriesNum);
                            disp(getReport(ME));
                        end
                        v_scale=PixelSize/str2double(SeriesInfo{71,2});
                        
                        for Ipart=1:size(PartInfoC{ii},1)
                            Itot=Itot+1;
                            info{Itot,5} = SeriesNum; %#ok<*AGROW>
                            
                            %                         SeriesInd=find(AllSeriesInd==info{Itot,5});
                            %                         info{Itot,6} = SeriesTime(SeriesInd);
                            info{Itot,7} = str2double(SeriesInfo{76,2});%GainValue
                            info{Itot,8} = str2double(SeriesInfo{72,2});%EMgain
                            info{Itot,9} = str2double(SeriesInfo{71,2});%CicleTime (delta_t?)
                            info{Itot,10} = str2double(SeriesInfo{end,2});%NumberOfTimeStamp
                            a=PartInfoC{ii}{Ipart,5};
                            info{Itot,15} = a; %risultati fit
                            info{Itot,16} = TimeFromAd(strcmp(A.textdata,[FolderList(Ifolder).name,'.lif']) & SeriesPos==SeriesNum);
                            if ~calc_exC && isfield(extrC{ii}(Ipart),'orientation') && ~isempty(extrC{ii}(Ipart).orientation)
                                orientation=extrC{ii}(Ipart).orientation;
                            else
                                orientation='x';
                                extrC{ii}(Ipart).orientation=orientation;
                            end
                            switch(orientation(1))
                                case 'l'
                                    ornv=[-1,0];
                                case 'u'
                                    ornv=[0,-1];
                                case 'd'
                                    ornv=[0,1];
                                otherwise
                                    ornv=[1,0];
                            end
                            extrC{ii}(Ipart).ornv=ornv;
                            if calc_exC || calcexc
                                a=a(repmat(a(:,3)<4,[1 3]));
                                a=reshape(a,[numel(a)/3 3]); %solo "raggi" < 4?
                                info{Itot,2}=mean(a(:,1).*a(:,3).^2); %intensità media?
                                Y=PartInfoC{ii}{Ipart,2};
                                X=PartInfoC{ii}{Ipart,1};
                                T=PartInfoC{ii}{Ipart,3};
                                %from Fase2newLuin start
                                x_av=movingmean(X,av_win);
                                y_av=movingmean(Y,av_win);
                                t_av=movingmean(T,av_win);

                                switch(orientation(1))
                                    case 'l'
                                        orn=-sign(vx_av);
                                    case 'u'
                                        orn=-sign(vy_av);
                                    case 'd'
                                        orn=sign(vy_av);
                                    otherwise
                                        orn=sign(vx_av);
                                end
                                vx_av=diff(x_av)./diff(t_av);
                                vy_av=diff(y_av)./diff(t_av);
                                v_av=sqrt(vx_av.^2+vy_av.^2);
                                v_av=v_av.*orn;
                                if ~tip || orientation(1)=='x'
                                    vv_av=vx_av;
                                else
                                    vv_av=v_av;
                                end
                                [S,starts,nsub]=StopFwRv(vv_av,vmin*info{Itot,9}/PixelSize,tc,t_av);
                                extrC{ii}(Ipart).SelId(Ipart)=Itot;
                                extrC{ii}(Ipart).SelIdInd=1:length(T);
                                extrC{ii}(Ipart).x=X;
                                extrC{ii}(Ipart).y=Y;
                                extrC{ii}(Ipart).t=T;
                                
                                [extrC{ii}(Ipart).xmin(1), extrC{ii}(Ipart).xmin(2)]=min(extrC{ii}(Ipart).x);
                                %                            extrC{ii}(Ipart).min(2)=extrC{ii}(Ipart).t(b);
                                [extrC{ii}(Ipart).xmax(1), extrC{ii}(Ipart).xmax(2)]=max(extrC{ii}(Ipart).x);
                                %                            extrC{ii}(Ipart).max(2)=ind(b);
                                extrC{ii}(Ipart).vx=diff(extrC{ii}(Ipart).x)./diff(extrC{ii}(Ipart).t);
                                extrC{ii}(Ipart).vy=diff(extrC{ii}(Ipart).y)./diff(extrC{ii}(Ipart).t);
                                extrC{ii}(Ipart).t_av=t_av;
                                extrC{ii}(Ipart).x_av=x_av;
                                extrC{ii}(Ipart).y_av=y_av;
                                extrC{ii}(Ipart).vx_av=vx_av;
                                extrC{ii}(Ipart).vy_av=vy_av;
                                extrC{ii}(Ipart).v_av=v_av;
                                extrC{ii}(Ipart).S=S;
                                extrC{ii}(Ipart).starts=starts;
                                extrC{ii}(Ipart).nsub=nsub;
                                %from Fase2newLuin end
                            else
                                Y=extrC{ii}(Ipart).y;
                                X=extrC{ii}(Ipart).x;
                                T=extrC{ii}(Ipart).t;
                                
                            end
                            info{Itot,17} = T;
                            info{Itot,11} = Y;
                            info{Itot,12} = X;
                            X=SmoothGauss(X,2);
                            Y=SmoothGauss(Y,2);
                            T=SmoothGauss(T,2);
                            dx=X(2:end)-X(1:end-1);
                            dy=Y(2:end)-Y(1:end-1);
                            [Th,Rho]=cart2pol(dx,dy);
                            info{Itot,13} = Th; %angolo di spostamento
                            info{Itot,14} = Rho; %spostamento
                            
                            D=sqrt((X(1)-X).^2+(Y(1)-Y).^2);
                            info{Itot,1}=sign(X(1)-X(end))*(D(:))*0.228;%spostamento medio con segno
                            D=SmoothGauss(D,1);
                            %V=(D(2:end)-D(1:end-1))/(T(2:end)-T(1:end-1))/str2double(SeriesInfo{71,2}); %NON USATA! velocità "istantanea"?; 1 ERRORE rimasto: considera differenze fra distanze dalla partenza
                            info{Itot,3}=sign(X(1)-X(end))*(D(end))/(T(end)-T(1))/str2double(SeriesInfo{71,2})*0.228; %velocità media con segno?
                            D=sqrt(X.^2+Y.^2);
                            D=(D(1)-D(:))*0.228; %??? differenza fra moduli di vettore invece che fra vettori?!?
                            D=SmoothGauss(D,2);
                            V=(D(2:end)-D(1:end-1))/str2double(SeriesInfo{71,2}); %nel calcolare la differenza, il D(1) se ne va, corrisponde alla velocità con cui la particella si avvicina all'origine.
                            %                 V1=[V1; V];
                            if sum(V>0.5)>1 %any(V>0.5)???
                                info{Itot,4}=mean(V(V>0.5)); %forse vecchio metodo per calcolare velocità media su pezzi non fermi.
                            end
                            
                            %inizio nuovo Luin
                            %                        xC{ii}=extrC{ii}(Ipart);
                            [xmin, ixmin]=min(extrC{ii}(Ipart).x_av);
                            [xmax, ixmax]=max(extrC{ii}(Ipart).x_av);
                            [ymin, iymin]=min(extrC{ii}(Ipart).y_av);
                            [ymax, iymax]=max(extrC{ii}(Ipart).y_av);
                            dxmax=(xmax-xmin)*PixelSize;
                            dymax=(ymax-ymin)*PixelSize;
                            
                            if ornv(1)
                                extrC{ii}(Ipart).vxmax=((xmax-xmin)/(extrC{ii}(Ipart).t_av(ixmax)-extrC{ii}(Ipart).t_av(ixmin)))*v_scale;
                                vzmax=ornv(1)*extrC{ii}(Ipart).vxmax;
                                vzperp=(v_scale*(extrC{ii}(Ipart).y_av(ixmax)-extrC{ii}(Ipart).y_av(ixmin))/(extrC{ii}(Ipart).t_av(ixmax)-extrC{ii}(Ipart).t_av(ixmin)));
                                extrC{ii}(Ipart).vymax=vzperp;
                                dmax=sign(ornv(1)*dxmax)*sqrt(dxmax^2+(extrC{ii}(Ipart).y_av(ixmax)-extrC{ii}(Ipart).y_av(ixmin))^2*PixelSize^2);
                            else
                                extrC{ii}(Ipart).vymax=((ymax-ymin)/(extrC{ii}(Ipart).t_av(iymax)-extrC{ii}(Ipart).t_av(iymin)))*v_scale;
                                vzmax=ornv(2)*extrC{ii}(Ipart).vymax;
                                vzperp=(v_scale*(extrC{ii}(Ipart).x_av(iymax)-extrC{ii}(Ipart).x_av(iymin))/(extrC{ii}(Ipart).t_av(iymax)-extrC{ii}(Ipart).t_av(iymin)));
                                extrC{ii}(Ipart).vxmax=vzperp;
                                dmax=sign(ornv(2)*dymax)*sqrt(dymax^2+(extrC{ii}(Ipart).x_av(iymax)-extrC{ii}(Ipart).x_av(iymin))^2*PixelSize^2);
                            end
                            extrC{ii}(Ipart).vmax=sign(vzmax)*sqrt(vzmax^2+vzperp^2);
                            
                            extrC{ii}(Ipart).vx_old=((extrC{ii}(Ipart).x_av(end)-extrC{ii}(Ipart).x_av(1))/(extrC{ii}(Ipart).t_av(end)-extrC{ii}(Ipart).t_av(1)))*v_scale;
                            extrC{ii}(Ipart).vy_old=((extrC{ii}(Ipart).y_av(end)-extrC{ii}(Ipart).y_av(1))/(extrC{ii}(Ipart).t_av(end)-extrC{ii}(Ipart).t_av(1)))*v_scale;
 
                            extrC{ii}(Ipart).v_old=sign(ornv*[extrC{ii}(Ipart).vx_old;extrC{ii}(Ipart).vy_old])...
                                *sqrt(extrC{ii}(Ipart).vx_old^2+extrC{ii}(Ipart).vy_old^2);
                            
                            extrC{ii}(Ipart).dxmax=dxmax;
                            extrC{ii}(Ipart).dymax=dymax;
                            extrC{ii}(Ipart).dmax=dmax;
                            extrC{ii}(Ipart).t_tot=(extrC{ii}(Ipart).t(end)-extrC{ii}(Ipart).t(1)+1)*str2double(SeriesInfo{71,2}); % t in frames
                            extrC{ii}(Ipart).n_tot=length(extrC{ii}(Ipart).t);
                            %dx=diff(extrC{ii}(Ipart).x_av(extrC{ii}(Ipart).starts))*PixelSize;
                            %dy=diff(extrC{ii}(Ipart).y_av(extrC{ii}(Ipart).starts))*PixelSize;
                            dt=diff(extrC{ii}(Ipart).t(extrC{ii}(Ipart).starts))*str2double(SeriesInfo{71,2});
                            dt(end)=dt(end)+str2double(SeriesInfo{71,2});
                            dn=diff(extrC{ii}(Ipart).starts);
                            extrC{ii}(Ipart).dt=dt; %sec
                            extrC{ii}(Ipart).dn=dn;
                            %                         extrC{ii}(Ipart).vxp=dx./dt;
                            %                         extrC{ii}(Ipart).vyp=dy./dt;
                            extrC{ii}(Ipart).vxp=zeros(size(dt));
                            extrC{ii}(Ipart).vyp=zeros(size(dt));
                            extrC{ii}(Ipart).vxp_s=zeros(size(dt));
                            extrC{ii}(Ipart).vyp_s=zeros(size(dt));
                            for Ic=1:length(dt) 
                                if dn(Ic)>3
                                    [a,b]=linfit(extrC{ii}(Ipart).t((extrC{ii}(Ipart).starts(Ic)+1):(extrC{ii}(Ipart).starts(Ic+1)-1)),extrC{ii}(Ipart).x((extrC{ii}(Ipart).starts(Ic)+1):(extrC{ii}(Ipart).starts(Ic+1)-1)));
                                    extrC{ii}(Ipart).vxp(Ic)=a(1)*v_scale;
                                    extrC{ii}(Ipart).vxp_s(Ic)=sqrt(b(1,1))*v_scale;
                                    [a,b]=linfit(extrC{ii}(Ipart).t((extrC{ii}(Ipart).starts(Ic)+1):(extrC{ii}(Ipart).starts(Ic+1)-1)),extrC{ii}(Ipart).y((extrC{ii}(Ipart).starts(Ic)+1):(extrC{ii}(Ipart).starts(Ic+1)-1)));
                                    extrC{ii}(Ipart).vyp(Ic)=a(1)*v_scale;
                                    extrC{ii}(Ipart).vyp_s(Ic)=sqrt(b(1,1))*v_scale;
                                else
                                    [a,b]=linfit(extrC{ii}(Ipart).t((extrC{ii}(Ipart).starts(Ic)):(extrC{ii}(Ipart).starts(Ic+1))),extrC{ii}(Ipart).x((extrC{ii}(Ipart).starts(Ic)):(extrC{ii}(Ipart).starts(Ic+1))));
                                    extrC{ii}(Ipart).vxp(Ic)=a(1)*v_scale;
                                    extrC{ii}(Ipart).vxp_s(Ic)=sqrt(b(1,1))*v_scale;
                                    [a,b]=linfit(extrC{ii}(Ipart).t((extrC{ii}(Ipart).starts(Ic)):(extrC{ii}(Ipart).starts(Ic+1))),extrC{ii}(Ipart).y((extrC{ii}(Ipart).starts(Ic)):(extrC{ii}(Ipart).starts(Ic+1))));
                                    extrC{ii}(Ipart).vyp(Ic)=a(1)*v_scale;
                                    extrC{ii}(Ipart).vyp_s(Ic)=sqrt(b(1,1))*v_scale;
                                end
                            end
                            [a,b]=linfit(extrC{ii}(Ipart).t,extrC{ii}(Ipart).x);
                            extrC{ii}(Ipart).vxmed=a(1)*v_scale;
                            extrC{ii}(Ipart).vxmed_s=sqrt(b(1,1))*v_scale;
                            [a,b]=linfit(extrC{ii}(Ipart).t,extrC{ii}(Ipart).y);
                            extrC{ii}(Ipart).vymed=a(1)*v_scale;
                            extrC{ii}(Ipart).vymed_s=sqrt(b(1,1))*v_scale;
                            extrC{ii}(Ipart).vmed=sign(ornv*[extrC{ii}(Ipart).vxmed;extrC{ii}(Ipart).vymed])*sqrt(extrC{ii}(Ipart).vxmed.^2+extrC{ii}(Ipart).vymed.^2);
                            extrC{ii}(Ipart).vmed_s=sqrt(extrC{ii}(Ipart).vxmed_s.^2.*extrC{ii}(Ipart).vxmed.^2+extrC{ii}(Ipart).vymed_s.^2.*extrC{ii}(Ipart).vymed.^2)./abs(extrC{ii}(Ipart).vmed);
                            
                            extrC{ii}(Ipart).vp=sqrt(extrC{ii}(Ipart).vxp.^2+extrC{ii}(Ipart).vyp.^2);
                            extrC{ii}(Ipart).vp_dir=sign(ornv*[extrC{ii}(Ipart).vxp;extrC{ii}(Ipart).vyp]).*extrC{ii}(Ipart).vp;
                            extrC{ii}(Ipart).vp_s=sqrt(extrC{ii}(Ipart).vxp_s.^2.*extrC{ii}(Ipart).vxp.^2+extrC{ii}(Ipart).vyp_s.^2.*extrC{ii}(Ipart).vyp.^2)./extrC{ii}(Ipart).vp;
                            [extrC{ii}(Ipart).v_avsn,extrC{ii}(Ipart).v_avsn_s2m,extrC{ii}(Ipart).v_avsn_s2s]=wmean(extrC{ii}(Ipart).vp(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0),1./(extrC{ii}(Ipart).vp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0)).^2,[],1);
                            [extrC{ii}(Ipart).v_avdx,extrC{ii}(Ipart).v_avdx_s2m,extrC{ii}(Ipart).v_avdx_s2s]=wmean(extrC{ii}(Ipart).vp(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0),1./(extrC{ii}(Ipart).vp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0)).^2,[],1);
                            [extrC{ii}(Ipart).v_aval,extrC{ii}(Ipart).v_aval_s2m,extrC{ii}(Ipart).v_aval_s2s]=wmean(extrC{ii}(Ipart).vp(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0),1./(extrC{ii}(Ipart).vp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0)).^2,[],1);
                            %IMPORTANT: next ones are still called 'vx' for
                            %code compability, but they are actually the
                            %total velocity oriented with respect to the
                            %given orientation. vx_avsn is still positive
                            %and not negative, but vx_aval now is with sign.
                            [extrC{ii}(Ipart).vx_avsn,extrC{ii}(Ipart).vx_avsn_s2m,extrC{ii}(Ipart).vx_avsn_s2s]=wmean(extrC{ii}(Ipart).vp_dir(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0),1./(extrC{ii}(Ipart).vp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0)).^2,[],1); 
                            extrC{ii}(Ipart).vx_avsn=-extrC{ii}(Ipart).vx_avsn;
                            [extrC{ii}(Ipart).vx_avdx,extrC{ii}(Ipart).vx_avdx_s2m,extrC{ii}(Ipart).vx_avdx_s2s]=wmean(extrC{ii}(Ipart).vp_dir(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0),1./(extrC{ii}(Ipart).vp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0)).^2,[],1);
                            [extrC{ii}(Ipart).vx_aval,extrC{ii}(Ipart).vx_aval_s2m,extrC{ii}(Ipart).vx_aval_s2s]=wmean(extrC{ii}(Ipart).vp_dir(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0),1./(extrC{ii}(Ipart).vp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0)).^2,[],1);
                            % [extrC{ii}(Ipart).vx_avsn,extrC{ii}(Ipart).vx_avsn_s2m,extrC{ii}(Ipart).vx_avsn_s2s]=wmean(extrC{ii}(Ipart).vxp(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0),1./(extrC{ii}(Ipart).vxp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0)).^2,[],1); 
                            % [extrC{ii}(Ipart).vx_avdx,extrC{ii}(Ipart).vx_avdx_s2m,extrC{ii}(Ipart).vx_avdx_s2s]=wmean(extrC{ii}(Ipart).vxp(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0),1./(extrC{ii}(Ipart).vxp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0)).^2,[],1);
                            % [extrC{ii}(Ipart).vx_aval,extrC{ii}(Ipart).vx_aval_s2m,extrC{ii}(Ipart).vx_aval_s2s]=wmean(extrC{ii}(Ipart).vxp(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0),1./(extrC{ii}(Ipart).vxp_s(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0)).^2,[],1);
                            extrC{ii}(Ipart).dtsn=sum(dt(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0));
                            extrC{ii}(Ipart).dtdx=sum(dt(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0));
                            extrC{ii}(Ipart).dtal=sum(dt(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0));
                            extrC{ii}(Ipart).dnsn=sum(dn(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))<0));
                            extrC{ii}(Ipart).dndx=sum(dn(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))>0));
                            extrC{ii}(Ipart).dnal=sum(dn(extrC{ii}(Ipart).S(extrC{ii}(Ipart).starts(1:end-1))==0));
                            extrC{ii}(Ipart).trname=[FolderList(Ifolder).name ' ' FileName(1:end-4) ' ' int2str(SeriesNum) ' C' int2str(ii) ' ' int2str(Ipart)];
                            extrC{ii}(Ipart).t_from_add=info{Itot,16};
                        end
                        extr=catstruct(extr,extrC{ii});
                    end
                end
                
                if calcexc, FileName=['new',FileName]; end
                extrC1=extrC{1}; %#ok<NASGU>
                extrC2=extrC{2}; %#ok<NASGU>
                PartInfoC1=PartInfoC{1}; %#ok<NASGU>
                PartInfoC2=PartInfoC{2}; %#ok<NASGU>
                save(FileName,'PartInfoC2','PartInfoC1','extrC1','extrC2');
            end
        end
        cd('..')
    end
end

if calcexc, appx=['new_05',appx]; end %?
if tip, appx=['tip',appx]; end

save([InfoName appx '.mat'],'info','extr');
warning(stat);
disp(['Fase3 end ' datestr(now)]);