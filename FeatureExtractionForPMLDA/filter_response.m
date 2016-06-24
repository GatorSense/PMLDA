function [out]=filter_response(I,Seg)
%This function is to compute the filter response for sonar image
%Input:
%   I  - input image of size M-by-N
%   Seg - Segmentation map (M-by-N)
%Output:
%   out - filter response 


[m n]=size(I);
dXImg=0:0.0125:(n-1)*0.0125;
dXImg=repmat(dXImg,m,1);

Labels=unique(Seg);
Response=[];
for i=1:length(Labels)
i
        [row col]=find(Seg==i);
        idx=sub2ind(size(Seg),row,col);
        clip=I.*(Seg==i);
        
        [freqBest, ABest, bBest,D0] = metropolisSandRipple(clip,dXImg, 2,Seg==i);
        filter_len=round(40/freqBest);
        LL=round(filter_len/3);
        neg_1 = -1*ones(11,LL);
        pos_1 = ones(11,LL);
        zero_p = zeros(11,LL);
        repeat_pattern = [neg_1, zero_p, pos_1, zero_p, neg_1, zero_p, pos_1];
        out1=conv2(I,repeat_pattern,'same'); 
        Response=[Response abs(out1(:))];

end


out = zeros(size(Seg));
for i=1:length(unique(Labels))  
    res=Response(:,i);
    res=reshape(res,size(I,1),size(I,2));
    res=res.*(Seg==i);
    out = out + res;
end
