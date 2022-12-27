function [S,starts,nsub]=StopFwRv(v,vc,tc,t)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
if nargin<1 || isempty(v)
    S=NaN;
    starts=NaN;
    nsub=NaN;
    return
end
sv=size(v);
lv=length(v);
v=v(:)';

if nargin<4 || isempty(t) || (length(t)~=lv && length(t)~=lv+1)
    t=1:lv+1;
elseif length(t)==lv
    t(end+1)=t(end)+(t(end)-t(1))/lv;
end
t=t(:)';
dt=diff(t);

if nargin<3 || isempty(tc)
    tc=1;
end
if tc<0, tc=-tc; end

if nargin<2 || isempty(vc)
    vc=1;
end
if vc<0, vc=-vc; end

if lv<2*tc
    vav=(v*dt')/sum(dt);
    S(sv)=(vav>vc)-(vav<-vc);
    nsub=1;
    starts=[1,lv+1];
    return
end

S=(v>vc)-(v<-vc);
[starts,nsub]=STALL_divide_sub(S);
dstart=diff(starts);
stw=find(dstart<tc);
if ~isempty(stw)
    dstw=diff(stw);
    if any(dstw==1)
        [dstwstarts, nseg, isin]=STALL_divide_sub(dstw==1);
        for II=1:nseg
            if isin<0
                range=starts(stw(dstwstarts(II))):(starts(stw(dstwstarts(II+1))+1)-1);
                %                 vav=(v(range)*dt(range)')/sum(dt(range)); %mean(v(range));
                vav=mov_mean_all(v(range),tc,dt(range));
                SS=(vav>vc)-(vav<-vc);
                Sp1=find(SS==1);
                Sm1=find(SS==-1);
                Sp=Sp1;
                Sm=Sm1;
                for ttc=1:(tc-1)
                    Sp1=unique([Sp1,Sp+ttc]);
                    Sm1=unique([Sm1,Sm+ttc]);
                end
                Sp1=Sp1(Sp1<=length(range));
                Sm1=Sm1(Sm1<=length(range));
                Spm1=intersect(Sp1,Sm1);
                S(range)=0;
                S(range(Sp1))=1;
                S(range(Sm1))=-1;
                %Spm1bigger0=v(range(Spm1))>0;
                %S(range(Spm1(Spm1bigger0)))=1;
                S(range(Spm1))=(v(range(Spm1))>vc)-(v(range(Spm1))<-vc);
                startsl=STALL_divide_sub(S(range));
                dstartl=diff(startsl);
                stwl=find(dstartl>tc);
                if ~isempty(stwl) %let's try to see if I can reduce the size of "moving" parts (only possible if Sp/m1 had subsequents indexes)
                    for errc=1:length(stwl)
                        rangel=startsl(stwl(errc)):(startsl(stwl(errc)+1)-1);
                        Sin=S(range(rangel(1)));
                        if Sin
                            lr=length(rangel);
                            if range(rangel(1))==1 || S(range(rangel(1))-1)
                                if range(rangel(end))~=lv && ~S(range(rangel(end))+1)
                                    for ccl=0:(lr-tc-1)
                                        if Sin*v(range(rangel(end-ccl)))<vc
                                            S(range(rangel(end-ccl)))=0;
                                        else
                                            break
                                        end
                                    end
                                end
                            else
                                if range(rangel(end))~=lv && ~S(range(rangel(end))+1)
                                    for ccl=1:(lr-tc)
                                        [vminl,begend]=min([Sin*v(range(rangel(1))),Sin*v(range(rangel(end)))]);
                                        if vminl<vc
                                            if begend==1
                                                S(range(rangel(1)))=0;
                                                rangel(1)=[];
                                            else
                                                S(range(rangel(end)))=0;
                                                rangel(end)=[];
                                            end
                                        else
                                            break
                                        end
                                    end
                                else
                                    for ccl=1:(lr-tc)
                                        if Sin*v(range(rangel(ccl)))<vc
                                            S(range(rangel(ccl)))=0;
                                        else
                                            break
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            isin=-isin;
        end
        [starts,nsub]=STALL_divide_sub(S);
        dstart=diff(starts);
        stw=find(dstart<tc);
    end
    
    if ~isempty(stw)
        dstw=diff(stw);
        if any(dstw==1)
            [dstwstarts, nseg, isin]=STALL_divide_sub(dstw==1);
            for II=1:nseg
                if isin<0
                    range=starts(stw(dstwstarts(II))):(starts(stw(dstwstarts(II+1))+1)-1);
                    vav=(v(range)*dt(range)')/sum(dt(range)); %mean(v(range));
                    if vav>vc, S(range)=1;
                    elseif vav<-vc, S(range)=-1;
                    else S(range)=0;
                    end
                end
                isin=-isin;
            end
        end
        [starts,nsub]=STALL_divide_sub(S);
        dstart=diff(starts);
        stw=find(dstart<tc);
    end
    
    if ~isempty(stw)
        for errc=1:length(stw)
            range=starts(stw(errc)):(starts(stw(errc)+1)-1);
            lr=length(range);
            if stw(errc)==nsub
                [starts,nsub]=STALL_divide_sub(S);
                dstart=diff(starts);
                stw=find(dstart<tc);
                if dstart(end)<tc
                    if dstart(end-1)>2*tc-lr
                        rangeext=(min(range)-(tc-lr)):range(end);
                    else
                        rangeext=starts(end-2):range(end);
                    end
                    vav=(v(rangeext)*dt(rangeext)')/sum(dt(rangeext));
                    S(rangeext)=(vav>vc)-(vav<-vc);
                    %S(starts(stw(errc)):(starts(stw(errc)+1)-1))=S(starts(stw(errc)-1));
                end
            else
                if stw(errc)==1 || S(starts(stw(errc)+1))==S(starts(stw(errc)-1))
                    if abs(S(starts(stw(errc)+1))-S(starts(stw(errc))))<=1
                        S(range)=S(starts(stw(errc)+1));
                    else
                        Sout=S(starts(stw(errc)+1));
                        if stw(errc)==1
                            minrangeext=1;
                        else
                            minrangeext=max(starts(stw(errc)-1)+tc,1);
                        end
                        if stw(errc)+2>length(starts)
                            maxrangeext=length(v)-tc;
                        else
                            maxrangeext=starts(stw(errc)+2)-tc-1;
                        end
                        rangeext=max(min(range)-tc+lr,minrangeext):min(max(range)+tc-lr,maxrangeext);
                        lre=length(rangeext);
                        if lre<tc
                            if stw(errc)+2>length(starts)
                                maxdx=length(v);
                            else
                                maxdx=starts(stw(errc)+2)-1;
                            end
                            rangedx=starts(stw(errc)):maxdx;
                            vavdx=(v(rangedx)*dt(rangedx)')/sum(dt(rangedx));
                            if stw(errc)<=1
                                S(rangedx)=(vavdx>vc)-(vavdx<-vc);
                            else
                                minsn=starts(stw(errc)-1);
                                rangesndx=minsn:maxdx;
                                vav=(v(rangesndx)*dt(rangesndx)')/sum(dt(rangesndx));
                                if Sout*vav>vc
                                    S(range)=Sout;
                                else
                                    rangesn=minsn:starts(stw(errc)+1)-1;
                                    vavsn=(v(rangesn)*dt(rangesn)')/sum(dt(rangesn));
                                    [~,sndx]=min([Sout*vavsn,Sout*vavdx]);
                                    switch sndx
                                        case 1
                                            S(rangesn)=-Sout*(abs(vavsn)>vc);
                                        case 2
                                            S(rangedx)=-Sout*(abs(vavdx)>vc);
                                    end
                                end
                            end
                        else
                            vav=mov_mean_all(v(rangeext),tc,dt(rangeext));
                            if any(Sout*vav<vc)
                                [vav,rangeextst]=min(Sout*vav);
                                S(rangeext(rangeextst:min(rangeextst+tc-1,lre)))=-Sout*(abs(vav)>vc);
                            else
                                S(range)=Sout;
                            end
                        end
                    end
                else
                    if (v(range)*dt(range)')/sum(dt(range))>0
                        S(range)=max(S(starts(stw(errc)+1)),S(starts(stw(errc))-1));
                    else
                        S(range)=min(S(starts(stw(errc)+1)),S(starts(stw(errc))-1));
                    end
                end
            end
        end
    end
    
end
S=reshape(S,sv);
if nargout>1
    [starts,nsub]=STALL_divide_sub(S);
end





