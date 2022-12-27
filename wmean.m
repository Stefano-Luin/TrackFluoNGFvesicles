function [y,sigma2_m,sigma2_s,N,sigma2_mi,sigma2_si,y_si,W_tot] =...
    wmean(x,w,dim,chil1,sigma, s_type, w_type)
%[y,sigma2_m,sigma2_s,N,sigma2_mi,sigma2_si,y_si] =...
%    wmean(x,w,dim,chil1,sigma, s_type, w_type)
%   WMEAN   Weighted Average or mean value.
%   For vectors, WMEAN(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, WMEAN(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   WMEAN(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   WMEAN(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       wmean(x,w)
%
%   Luin: added calculation of variances:
%   sigma2_m: SE^2
%   sigma2_s: SD^2
%   N: number of acceptable (non 0 weight, no NaNs) "entries")
%   sigma2_mi,sigma2_si,y_si: results calculated based on equations from
%   [Gatz, D. F., and L. Smith. 1995 Atmospheric Environment 29:1185-1193]
%   (weighted mean errors) or [Rukhin 2007, Statistics & Probability
%   Letters 77:853–861] (confidence intervals with t=1 if sigma are given)
%
%   Other input:
%   chil1:  1:  check if sigma2_m might be underestimated (if calculated
%               reduced chisquare is less than 1), considering weighting
%               based on uncertainties (only if explicit). 
%           0:  standard if no w specified. SE=0 if it cannot be
%               calculated.
%   sigma:  uncertainties in the data; if present, w is considered a
%           number (default) or frequency of observations, depending on
%           w_type.
%   s_type: 0:  sigma are not considered.
%           1:  sigma are standard errors (default if present)
%           2:  sigma are standard deviations
%           3:  weighting becomes w(normalized if w_type 2, else 1)./sigma.^2,
%               also in sigma2_mi, sigma2_si (sigma SE by default)
%   w_type: 1:  consider w as number of observations, already included (or
%               to be included, depending on s_type=1 or 2) in the
%               calculation of standard errors in the data.
%           2:  consider w as propto frequency of observations, not at
%               all as an uncertainty^(-2).
%           3:  consider w as SE^-2
%           4:  consider w as SD^-2
%           5:  consider w as number of observations, considering them as
%               number of repetition of measurements giving exactly the
%               same results.
%           0:  unknow meaning of w. Treated as 2 or 3 depending on the
%               other specifications and on the output.
%             
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)

%% wmean
if nargin<1 || isempty(x)
    warning('Luin:wmean','Not enough input arguments.');
    res=num2cell(nan(1,7));
    [y,sigma2_m,sigma2_s,N,sigma2_mi,sigma2_si,y_si]=res{:};
    return;
end
if isempty(x)
    y=NaN;
    sigma2_m=Inf;
    sigma2_s=Inf;
    N=0;
    sigma2_mi=Inf;
    sigma2_si=Inf;
    y_si=NaN;
    return;
end
if nargin<7 || isempty(w_type)
    w_type=0;
end

if nargin<4 || isempty(chil1)
    chil1=1;
end

if nargin<3||isempty(dim)
    % Determine which dimension SUM will use
    dim=[];
    if exist('w','var') && ~isempty(w)
        dim = find(size(w)~=1, 1);
    end
    if isempty(dim)
        dim = find(size(x)~=1, 1);
    end
    if isempty(dim), dim = 1; end
end

sizx=size(x);
if nargin<2 || isempty(w)
    w=ones(sizx);
    if nargin<2, chil1=0; end
    w_type = 5;
else
    sizw=size(w);
    % Check that dimensions of X match those of W.
    if(~isequal(sizx, sizw))
        if sizx(dim)==sizw(dim) && numel(w)==length(w)
            sizx(dim)=1;
            w=repmat(w,sizx);
        else
            error('Inputs x and w must be the same size.');
        end
    end
end
clear sizx sizw;
if nargin>=5 && ~isempty(sigma)
    % Check that dimensions of X match those of sigma.
    if(~isequal(size(x), size(sigma)))
        error('Inputs x and sigma must be the same size.');
    end
    if nargin<6 || isempty(s_type) || all(s_type ~= [0,1,2,3])
        s_type=1;
    end
    if all(w_type~=[1,2,5,0])
        warning('Luin:wmean','wrong weight definition if sigma is present; considering w_type=2');
        pause(2);
        w_type=2; % check...
    end
else
    if nargin<6 || isempty(s_type)
        s_type=0; % check!
    end
end
    
% Check that all of W are non-negative.
if (any(w(:)<0))
    warning('Luin:wmean','All weights, W, must be non-negative.');
    w(w<0)=0;
    pause(2);
end

isNaN=isnan(w)|isnan(x);
if s_type
    isNaN=isNaN|isnan(sigma);
    if w_type==1
        isNaN=isNaN|(w<=1);
    end
    sigma(isNaN)= Inf;
end
w(isNaN)=0;
x(isNaN)=0;
if nargin<5 || isempty(sigma)
    switch w_type
        case 3
            sigma=1./sqrt(w);
            s_type=1;
        case 4
            sigma=1./sqrt(w);
            s_type=2;
        case 0
            sigma=ones(size(x));
        otherwise
            sigma=zeros(size(x));
    end
end

if s_type==2
    if w_type == 1
        sigma=sigma./sqrt(w);
        sigma(isNaN)=Inf;
    else
        warning('Luin:wmean','wmean: error in specifications; considering sigma as SE instead of SD');
        s_type = 1;
        pause(2);
    end
end

% Check that there is at least one non-zero weight.
%if (all(w(:)==0))
%    error('At least one weight must be non-zero.');
%end


% N = size(x,dim);
N = sum(w>0,dim);
V1 = sum(w,dim);
W_tot=V1;
ni = ones(size(w)).*(w>0); %ni: number of observations in the population i
Ni = ni;            %Ni: number of repetition of observation i
fi = ni;            %fi: propto fraction if w_type==2
Ntot=N;
switch w_type
    case 1
        Ntot = N;
        ni = w;
    case 5
        Ntot = V1;
        Ni = w;
    case 2
        fi = bsxfun(@times,bsxfun(@rdivide,w,sum(w,dim)),N);
end

Nfull=sum(ni.*Ni,dim);

if s_type && s_type ~= 3
    switch w_type
        case {1,5}
            w=1./sigma.^2;
        case 0
            w=w./sigma.^2;
        case 2
            w=fi./sigma.^2;
    end
end
if s_type==3
    w=fi./sigma.^2;
    if w_type~=2
        warning('Luin:wmean','wmean: relayable results if s_type=3 only if w_type=2');
    end
end

V1 = sum(Ni.*w,dim);
w_i = bsxfun(@rdivide,w,V1);

if (any(N<2))
    warning('LUIN:wmean1','WMEAN: 1 or 0 non-zero weight for at least one entry');
end

y = sum(Ni.*w_i.*x,dim);
if nargout>6, y_si=y;end

%V1=sum(w,dim);
Chi2=sum(Ni.*w.*(bsxfun(@minus,x,y).^2),dim);
sigma2_m = (1./V1).*Chi2./(Ntot-1);
sigma2_s = V1./(V1.^2-sum(Ni.*w.^2,dim)).*Chi2;
sigma2_sem = 1./sum(fi.*Ni./sigma.^2,dim);

if nargout > 4
    Ntemp=[];
    sigma2_mi=sigma2_sem;
    if nargout > 5, sigma2_si=sigma2_sem.*Nfull; end
    if s_type==3
        fi=ones(size(fi)); %already included in w
        Ntemp = N;
        N = Ntot;
        s_type=0;
    end
    if s_type
        if w_type == 1
            sigma2_mi = sigma2_sem.* ... %1/sum(1./sigma.^2,dim)
                sum(arrayfun(@(c,z) hypergeomloc([1,2],c,z),(ni+1)/2,1-w_i).*w_i.*fi.*Ni,dim);
        else
            ptot=sum(Ni,dim);
            yDL=(sum(Ni.*fi.*(bsxfun(@minus,x,y).^2)./sigma.^2,dim)-ptot+1)./...
                (1./sigma2_sem-sum(Ni.*fi./sigma.^4,dim).*sigma2_sem);
            yDL=max(yDL,0);
            w_i=bsxfun(@plus,sigma,yDL).^(-2);
            w_i=bsxfun(@rdivide,w_i,sum(Ni.*fi.*w_i,dim));
            sigma2_mi = sum(Ni.*fi.*w_i.*(bsxfun(@minus,x,y).^2),dim)./...
                ((ptot-1).*prod(bsxfun(@(wi,p) (p.*wi).^(1./(p-1)),w_i,ptot).^(Ni.*fi),dim));
            if nargout>6
                y_si=sum(Ni.*fi.*w_i.*x,dim);
            end
        end
        if(chil1)
            sigma2_mi(sigma2_mi<sigma2_sem)=sigma2_sem(sigma2_mi<sigma2_sem);
        end
        if nargout > 5
            sigma2_si = sigma2_mi.*Nfull; %underestimated if w_type~=1
        end
    elseif any(w_type==[1,2,5,0])
        Pav=sum(w,dim)./N; %Pav=sum(w.^2,dim)./sum(w,dim);
        sigma2_mi = N./((N-1).*(sum(w,dim).^2)) .* (...
            sum((bsxfun(@minus,w.*x,Pav.*y).^2),dim) ...
            -2.*y.*sum(bsxfun(@minus,w,Pav).*bsxfun(@minus,w.*x,Pav.*y),dim) ...
            + y.^2.*sum(bsxfun(@minus,w,Pav).^2,dim));
        sigma2_si=Chi2./V1.*N./(N-1);
    end
    if ~isempty(Ntemp)
        N=Ntemp;
        s_type=3;
        %fi no more used
    end
end
if (chil1==1 && nargin>3)
    sigma2_m(Chi2<(Ntot-1))=1./V1(Chi2<(Ntot-1));
%   sigma2_m(Chi2<1)=1./V1(Chi2<1);
end

if any(N==0)
    av=sum(x,dim)./sum(~isNaN,dim);
    y(N==0) = av(N==0);
    sigma2_m(N==0) = Inf;
    sigma2_s(N==0) = Inf;
end

if nargin == 1
    sigma2_m(N==1) = Inf;
    sigma2_s(N==1) = Inf;
else
    if chil1==0
        if s_type
            sn0=sum(sigma.*(w~=0),dim);
            sigma2_m(N==1) = sn0(N==1).^2;
        else
            sigma2_m(N==1) = Inf; %check...
        end
    else
        sigma2_m(N==1) = 1./V1(N==1);
    end
    sigma2_s(N==1) = 1./V1(N==1);
end

if nargout>4
    if any(N==0)
%        av=sum(x,dim)./sum(~isNaN,dim);
        y_si(N==0) = av(N==0);
        sigma2_mi(N==0) = Inf;
        sigma2_si(N==0) = Inf;
    end

    if nargin == 1
        sigma2_mi(N==1) = Inf;
        sigma2_si(N==1) = Inf;
    else
        if chil1==0
            if s_type
 %               sn0=sum(sigma.*(w~=0),dim);
                sigma2_mi(N==1) = sn0(N==1).^2;
            else
                sigma2_mi(N==1) = Inf; %check...
            end
        else
            sigma2_mi(N==1) = 1./V1(N==1);
        end
        if s_type || ~w_type
            sigma2_si(N==1) = 1./V1(N==1);
        end
    end
end
end

% function d=hypergeomloc(a,b,c)
%         try
%             d=hypergeom(a,b,c);
%         catch ME
%             disp(getReport(ME,'extended'));
%             d=Inf;
%         end
% end
function d=hypergeomloc(a,b,c)
    d=zeros(size(b));
    cc=0;
    for e=b
        cc=cc+1;
        try
            d(cc)=hypergeom(a,b,c);
        catch ME
            disp(getReport(ME,'extended'));
            d(cc)=Inf;
        end
    end
end

%% test
% ll=[1,1,1; 2,2,2; 3,3,3] %#ok<UNRCH>
% w=[0.5,0.5,0;0.5,0,0;0.5,0,0]
% %w=[0.001,0.001,0;0.001,0,0;0.001,0,0]
% [y,sigma2_m,sigma2_s,N,smi,ssi,smw,ssw] = wmean_ongoing(ll,w,1,1)

