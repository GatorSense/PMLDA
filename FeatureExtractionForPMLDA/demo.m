%This code is to compute feature for sonar image and generate documents for
%the sonar image
%The variable SegData is a cell. Each cell represents a document
%SegData{i}.X - contains feature in document i
%SegData{i}.index - contains location of each pixel in document i

clc; clear all;close all
addpath frequencyHeightEstimation

load ImageHF-0-0_orig_seg
filter_size=21;           %filter size for mean and entropy
cut=filter_size-1;        %cut off the left and right boundary by "cut" size

[out]=filter_response(I,Seg);
out1=ordfilt2(out,49,true(7));  %smooth the filter output
h=fspecial('average',21);
out2=filter2(h,out1);
out=out2;
out=out(:,cut:end-cut);

entropy=entropyfilt(I,ones(filter_size)); %compute the entropy
entropy=entropy(:,cut:end-cut);

h=fspecial('average',filter_size); %compute the mean
Mean=filter2(h,I);
Mean=Mean(:,cut:end-cut);

I=I(:,cut:end-cut);
Seg=Seg(:,cut:end-cut);
ss=4;  %downsample the image and its features by ss
I=I(1:ss:end,1:ss:end);
Seg=Seg(1:ss:end,1:ss:end);
entropy=entropy(1:ss:end,1:ss:end);
Mean=Mean(1:ss:end,1:ss:end);
out=out(1:ss:end,1:ss:end);

figure(4);imagesc(I);colormap(gray);
figure(3);imagesc(entropy);title('entropy');hold on;
figure(2);imagesc(Mean);title('Mean');hold on;
figure(1);imagesc(out);title('Filter');hold on;

lc=Seg;
labels=unique(Seg);

for i=1:length(labels)  % split the image(feature)into documents according to the superpixel label
    i
    ll=labels(i);
    [row col]=find(lc==ll);
    index=sub2ind(size(I),row,col);
    f1=Mean(index);
    f2=entropy(index);
    f3=out(index);
    ff=[f1(:) f2(:) f3(:)];
    SegData{i}.X=ff;
    SegData{i}.index=index;
end



