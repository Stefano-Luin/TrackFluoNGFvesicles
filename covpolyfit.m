function cov=covpolyfit(S)
% Copyright (C) 2012 - 2022 Stefano Luin (s.luin@sns.it)
Rinv=inv(S.R);
cov=(Rinv*Rinv')*S.normr^2/S.df;
return
