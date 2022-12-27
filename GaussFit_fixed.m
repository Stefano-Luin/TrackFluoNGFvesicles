function [a, res,theor,ConfInt] = GaussFit_fixed(corr,pixelsize,X0,Y0,a0);
% Copyright (C) 2012 - 2022 Carmine di Rienzo - Stefano Luin (s.luin@sns.it)
h=1e-3;
set(gcbf,'pointer','watch');
% if ~exist('T0','var')
%     T0=0;
% end
% [X,Y] = meshgrid(-((size(corr,2)-1)/2)*pixelsize:pixelsize:((size(corr,2)-1)/2)*pixelsize,-((size(corr,1)-1)/2)*pixelsize:pixelsize:(size(corr,1)-1)/2*pixelsize);
[X,Y] = meshgrid(((1:size(corr,2))-X0)*pixelsize,((1:size(corr,1))-Y0)*pixelsize);
grid = [X Y];
% [~,R]=cart2pol(X,Y);
R = sqrt(X.^2+Y.^2);
%    
% [Y0, X0] = find(ismember(corr,max(max(corr))),size(corr,3));
% X0 = mod(X0,size(corr,2));
% 
% %EvilX0 and Y0 are where remainder from mod was zero -- these are set to
% %the "max" (ie size) of the corr
% EvilX0 = find(ismember(X0,0));
% X0(EvilX0) = size(corr,2);

    % Sets curve fit options, and sets lower bounds for amplitude and beam
    % radius to zero
%     lb = [0 0 -1 ];
%     ub = [];
    
    weights = ones(size(corr));
    
    % If there's whitenoise, 2 highest values (NB this might be more than
    % two points!) in corr func are set to zero, and given no weight in the fit
    

    

y0 = zeros(size(corr,3),1);
% g0 = max(reshape(corr(repmat(R~=0,[1 1 size(corr,3)])),[size(corr,1)*size(corr,2)-1 size(corr,3)]))';
g0=squeeze(corr(Y0,X0,:));
g0 = squeeze(g0);

% wguess = zeros(size(corr,3),1);
% for i=1:size(corr,3)
% [Wy, Wx] = find(ismember(abs((corr(:,:,i)/g0(i) - exp(-1))),min(min(abs(corr(:,:,i)/g0(i) - exp(-1))))));
% Wx = mod(Wx,size(corr,2));
% wguess(i) = mean(( (Wx - X0(i)).^2  + (Wy - Y0(i)).^2   ).^(1/2))*pixelsize;
% end

wguess = 0.4*ones(size(g0));

% % Converts from matrix index to LOCATION in pixelsize units
% for i=1:size(corr,3)
%   X0(i) = X(1,X0(i)); 
% end
% for i=1:size(corr,3)
%   Y0(i) = Y(Y0(i),1);  
% end

W=ones(size(corr,1),size(corr,2));
W(Y0,X0)=0;
W=W(:);
if exist('a0','var')
    initguess = a0;
else
    initguess = [g0 y0 wguess ];
end
    a=zeros(size(corr,3),3);
    ConfInt=zeros(size(corr,3),3,2);
    res=zeros(size(corr));
    theor=zeros(size(corr));
%     curvefitoptions = optimset('Display','off');
if size(corr,3)>1
    h = waitbar(0,'Fitting correlation functions...');
end
    for I=1:size(corr,3)
        ft = fittype( 'A2+A1.*exp(-(r.^2)./(A3.^(2)))', 'independent', {'r'},'dependent','z');
%         if T0==0 & I==1
%             W=ones(size(corr,1),size(corr,2));
%             W(Y0,:)=0;
%             W=W(:);
%         else
%             W=ones(size(corr,1),size(corr,2));
%             W(Y0,X0)=0;
%             W(Y0,:)=0;
%             corr=corr.*W;
%             W=W(:);
%         end
                opts = fitoptions(ft);
                opts.Display = 'Off';
%                 opts.MaxIter = 2000;
%                     opts.MaxFunEvals = 2000;
%                     opts.TolFun= 1.0000e-009;
%                 opts.Robust = 'On';
                opts.Lower = [0 0 0]; 
                opts.Upper = [40000 40000 5];
                if I==1
                    opts.StartPoint = initguess(I,:);
                else
                    opts.StartPoint = coeffvalues(fitresult);
                end
    %             opts.Upper = [Inf Inf Inf Inf +1.57];
                if I==1
                    W=ones(size(corr,1),size(corr,2));
%                     W(Y0,:)=0;
%                     W(Y0,X0)=0;
                    W=W(:);
                    opts.Weights=W;
                else
                    W=ones(size(corr,1),size(corr,2));
%                     W(Y0,X0)=0;
                    W=W(:);
                    opts.Weights=W;
                end
                 
            ACF=corr(:,:,I);
         [fitresult,gof,output] = fit( R(:), ACF(:), ft, opts );
         a(I,:)=coeffvalues(fitresult);
         if a(I,3)>5
             1
         end
         ConfInt(I,:,:)=confint(fitresult)';
         res(:,:,I)=reshape(output.residuals,[size(corr,1) size(corr,2)]);
         theor(:,:,I)=reshape(fitresult(R(:)),[size(corr,1) size(corr,2)]);
         if ishandle(h)
             
            waitbar(I/size(corr,3),h)
        end
         
    end
%     % Fits each corr func separately
%     if strcmp('2d',type)
%     for i=1:size(corr,3)
%         a0 = initguess(i,:);
%         a0xy(1:2) = initguess(i,1:2);
%         a0xy(3) = a0xy(2);
%         a0xy(4:6) = initguess(i,3:5);
%         [a(i,:),res(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit(@gauss2dwxy,a0xy,grid,corr(:,:,i).*weights(:,:,i),lb,ub,curvefitoptions,weights(:,:,i));
% 
%     end
%     end
% 
%     if strcmp('time',type)
%     for i=1:size(corr,3)
%         a0 = initguess(i,:);
%         [a(i,:),res(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit(@gauss2d,a0,grid,corr(:,:,i).*weights(:,:,i),lb,ub,curvefitoptions,weights(:,:,i));
%     end
%     end
if ishandle(h)
    close(h)
end
set(gcbf,'pointer','arrow');