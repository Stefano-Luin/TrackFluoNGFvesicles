function CrImg1=MovAv(CrImg,frame)
% Copyright (C) 2012 - 2022 Carmine di Rienzo (contact Stefano Luin s.luin@sns.it)


ImgM=mean(CrImg,3);
FrameNum=size(CrImg,3);



CrImg1=double(zeros(size(CrImg,1),size(CrImg,2),FrameNum-2*frame));
M=mean(CrImg(:,:,1:2*frame),3);
for i=frame:FrameNum-frame
    CrImg1(:,:,i-frame+1)=M;
    if i<FrameNum-frame
        M=M-double(CrImg(:,:,i-frame+1))/(2*frame)+double(CrImg(:,:,i+frame+1))/(2*frame);
    end
end
