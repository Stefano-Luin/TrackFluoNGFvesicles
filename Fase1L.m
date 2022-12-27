% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
clear
parDom=AskDomenica;
close all
FileList=dir;
for IFile=3:size(FileList,1)
% try
    if strcmp(FileList(IFile).name(end-2:end),'lif')
        FileName=FileList(IFile).name;
        mkdir(FileName(1:end-4))
        SizeOfChannel=128;
        info=ReadInfo(FileName,1,1);
%         FileName='Experiment005.lif';
        
        for I=1:info{6,2}
            if ~exist([FileName(1:end-4),filesep,'CH_Series_',num2str(I),'.mat'],'file')
                CrImg=bfopenLL_n(FileName,I,1,10000);
                Series.CrImg=CrImg{I};
                Series.ImgMeanT = squeeze(mean(mean(Series.CrImg,1),2));
                if any(Series.ImgMeanT==0)
                    Series.CrImg = Series.CrImg(:,:,1:(find(Series.ImgMeanT==0,1,'first')-1));
                end
                clear CrImg
                szImg=size(Series.CrImg);
                disp(size(Series.CrImg,3));
                if size(Series.CrImg,3)>50 && ~exist([FileName(1:end-4),filesep,'CH_Series_',num2str(I),'.mat'],'file')
                    disp(I);
                    Series.ImgMean=mean(Series.CrImg,3);
                    Channel.C1=[];
                    Channel.C2=[];
                    close all
                    imagesc(Series.ImgMean);axis image
                    N=input('Number of channel: ');
%%
                    if N==1 || N==2
                        h=imrect;
                        Pos=floor(wait(h)+1);
                        Pos(Pos<1)=1;
                        if Pos(2)+Pos(4)>szImg(1)
                            Pos(4)=szImg(1)-Pos(2);
                        end
                        if Pos(1)+Pos(3)>szImg(2)
                            Pos(3)=szImg(2)-Pos(1);
                        end
                        Channel.C1=Series.CrImg(Pos(2):Pos(2)+Pos(4),Pos(1):Pos(1)+Pos(3),:);
                        if N==2
                            h=imrect;
                            Pos=floor(wait(h)+1);
                            Pos(Pos<1)=1;
                            if Pos(2)+Pos(4)>szImg(1)
                                Pos(4)=szImg(1)-Pos(2);
                            end
                            if Pos(1)+Pos(3)>szImg(1)
                                Pos(3)=szImg(2)-Pos(1);
                            end
                            Channel.C2=Series.CrImg(Pos(2):Pos(2)+Pos(4),Pos(1):Pos(1)+Pos(3),:);
                        end
                    end
%%
                    save([FileName(1:end-4),filesep,'CH_Series_',num2str(I),'.mat'],'Channel')
                end
            end
        end
    end
% catch MEESnT1
%     disp(getReport(MEESnT1,'extended'));
%     save(['MEESnT1_',FileName(1:end-4),num2str(I),'.mat']);
%     continue;end
end
%%
for IFile=3:size(FileList,1)

    if strcmp(FileList(IFile).name(end-2:end),'lif')
        FileName=FileList(IFile).name;
        info=ReadInfo(FileName,1,1);
        for I=1:info{6,2}
                if exist([FileName(1:end-4),filesep,'CH_Series_',num2str(I),'.mat'],'file') && ~exist([FileName(1:end-4),filesep,'Series_',num2str(I),'.mat'],'file')
                    load([FileName(1:end-4),filesep,'CH_Series_',num2str(I),'.mat'])
                    Tracking.objsC1=[];
                    Tracking.objsC2=[];
                    Tracking.objs_linkC1=[];
                    Tracking.objs_linkC2=[];
                    'Tracking'
                    if ~isempty(Channel.C1)
                        ImgMeanT = squeeze(mean(mean(Channel.C1,1),2));
                    if sum(ImgMeanT==0)
                        Channel.C1 = Channel.C1(:,:,1:(find(ImgMeanT==0,1,'first')-1));
                        if ~isempty(Channel.C2)
                            Channel.C2 = Channel.C2(:,:,1:(find(ImgMeanT==0,1,'first')-1));
                        end
                    end

                        for j=1:size(Channel.C1,3)
                            [tmpobj] = fo5_rp(double(Channel.C1(:,:,j))-SmoothGauss(Channel.C1(:,:,j)...
                                ,parDom.GSmooth), parDom.processopt, parDom.processparam, parDom.thresh, parDom.fitstr);
                            tmpobj(5,:) = j;
                            Tracking.objsC1 = [Tracking.objsC1 tmpobj];
                        end
                            Tracking.objs_linkC1 = nnlink_rp(Tracking.objsC1,  parDom.step, parDom.memory);

                    end
                    if ~isempty(Channel.C2)
                            ImgMeanT = squeeze(mean(mean(Channel.C2,1),2));
                        if sum(ImgMeanT==0)
                            Channel.C2 = Channel.C2(:,:,1:(find(ImgMeanT==0,1,'first')-1));
                        end
                        for j=1:size(Channel.C2,3)
                            [tmpobj] = fo5_rp(double(Channel.C2(:,:,j))-SmoothGauss(Channel.C2(:,:,j)...
                                ,parDom.GSmooth), parDom.processopt, parDom.processparam, parDom.thresh, parDom.fitstr);
                            tmpobj(5,:) = j;
                            Tracking.objsC2 = [Tracking.objsC2 tmpobj];
                        end
                        Tracking.objs_linkC2 = nnlink_rp(Tracking.objsC2, parDom.step, parDom.memory);
                    end
                    save([FileName(1:end-4),filesep,'Series_',num2str(I),'.mat'],'Tracking')
                elseif exist([FileName(1:end-4),filesep,'Series_',num2str(I),'.mat'],'file')
                    load([FileName(1:end-4),filesep,'Series_',num2str(I),'.mat'])
                    
                    save([FileName(1:end-4),filesep,'Series_',num2str(I),'.mat'],'Tracking')
                end
                %%
         end
            
            close all
    end
end
