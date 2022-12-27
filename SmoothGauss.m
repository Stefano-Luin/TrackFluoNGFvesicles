function postConv=SmoothGauss(pCF,PSFSize)
% Copyright (C) 2012 - 2022 Carmine di Rienzo - Stefano Luin (s.luin@sns.it)

filtersize = ceil(PSFSize*6)+1;
postConv=zeros(size(pCF));

    
        for i1=1:size(pCF,3)
            postConv(:,:,i1) = convolveGaussian(squeeze(pCF(:,:,i1)),filtersize,PSFSize);
        end
