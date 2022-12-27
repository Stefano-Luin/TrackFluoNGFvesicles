function [starts,num_traj,isdiff]=STALL_divide_sub(b)
% STALL_divide_sub: script in STALL_divide;
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
sb=size(b,1);
if sb==1,b=b';end
starts=find([0;b]-[b;0]);
sizeb=size(b,1);
if isempty(starts)
    starts=[1;sizeb+1];
    isdiff=1;
else
    if starts(end)<=sizeb, starts(end+1)=sizeb+1; end
    if b(1),isdiff=-1;
    else
        isdiff=1;
        starts=[1;starts];
    end
end

num_traj=length(starts)-1;
if sb==1, starts=starts'; end
