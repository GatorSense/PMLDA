function [freqBest, ABest, bBest,D0_x,error] = metropolisSandRipple(clip, dXImg, ampSonar,Mask)

%Parameters
nmean         = 90;
nsig          = 50;
maxV          = 6.1e-3;
numIterations = 400;
dispRate      = 1000;
freqRange     = [3. 5];
ARange        = [0.1 0.2];

N = numel(clip);

%Sample Parameters
%For Frequency, Truncated Gaussian Mixture Model, centered on previous
%sample
mixF  = .7;
varBF = .1;
varSF = .01;

%For Amplitude, Truncated Gaussian Mixture Model, centered on previous
%sample
mixA  = .7;
varBA = .5;
varSA = .05;

%For Shift, Truncated Gaussian Mixture Model, centered on previous
%sample
mixS  = .7;
varBS = .25;
varSS = .01;

%%
%%%%%%%%%%%%
%Initialize%
%%%%%%%%%%%%

%Initialize Frequency, A and b
 % clip = medfilt2(clip, [9 9]);
[row col]=find(Mask>0);
mask=sub2ind(size(Mask),row,col);
clip_x=clip(mask);
clip_x=clip_x-mean(clip_x);
clip_x=clip_x/max(clip_x);


table=tabulate(row);
[maxCount,idx]=max(table(:,2));
most_row_idx=table(idx);
[coLL,~]=find(row==most_row_idx);
col1=col(coLL);
line=clip(most_row_idx,min(col1):max(col1));

[freq] = estimateFrequency(line, 1/(dXImg(1,2)-dXImg(1,1)));


 A = (1/2)*sqrt(max(clip_x(:)).^2+min(clip_x(:)).^2);

 A=0.05;
b = 0;


if(freq < freqRange(1))
    freq = freqRange(1);
elseif(freq > freqRange(2))
    freq = freqRange(2);
end

if(A < ARange(1))
    A = ARange(1);
elseif(A > ARange(2))
    A = ARange(2);
end

%%

D = simulateData(freq, A, b, size(clip,1), size(clip,2), dXImg(1,:), ampSonar);
D_x=D(mask);
D_x  = D_x - mean(D_x(:));
D_x = D_x / max(D_x(:));
objValue = (clip_x-D_x);
objValue = -(1/(N*2))*sum(sum(objValue.*objValue));

objValueBest = objValue;
freqBest = freq;
ABest  = A;
bBest = b;



for iter = 1:numIterations
    
    %disp(['iteration: ', num2str(iter), ' freqBest: ', num2str(freqBest), ' ABest: ', num2str(ABest), ' bBest: ', num2str(bBest)]);
    
    %%%%%
    %Sample Frequency
    %%%%%
%     display(['freqBest: ' num2str(freqBest)]);
%     display(['ABest: ' num2str(ABest)]);
   
    
    r = rand(1);
    if(r < mixF)
        freqSample = randn(1)*varSF + freq;
        while(freqSample < freqRange(1) || freqSample > freqRange(2))
            freqSample = randn(1)*varSF + freq;
        end
    else
        freqSample = randn(1)*varBF + freq;
        while(freqSample < freqRange(1) || freqSample > freqRange(2))
            freqSample = randn(1)*varBF + freq;
        end
    end
    
    %Evaluate Sample
    D = simulateData(freqSample, A, b, size(clip,1), size(clip,2), dXImg(1,:), ampSonar);
    D_x=D(mask);
    D_x  = D_x - mean(D_x(:));
    D_x = D_x / max(D_x(:));
    objValueS = (clip_x-D_x);
    objValueS = -(1/(N*2))*sum(sum(objValueS.*objValueS));
    
    diff = objValueS - objValue;
    P    = exp(diff);
    R    = rand;
    
    if(P >= R)
        freq = freqSample;
        objValue = objValueS;
  
        if(objValue > objValueBest)
            objValueBest = objValue;
            ABest = A;
            bBest = b;
            freqBest = freq;
        end
    end
    
    %%%%%
    %Sample A
    %%%%%
    r = rand(1);
    if(r < mixA)
        ASample = randn(1)*varSA + A;
        while(ASample < ARange(1) || ASample > ARange(2))
            ASample = randn(1)*varSA + A;
        end
    else
        ASample = randn(1)*varBA + A;
        while(ASample < ARange(1) || ASample > ARange(2))
            ASample = randn(1)*varBA + A;
        end
    end

    %Evaluate Sample
    D = simulateData(freq, ASample, b, size(clip,1), size(clip,2), dXImg(1,:), ampSonar);
    D_x=D(mask);
    D_x  = D_x - mean(D_x(:));
    D_x = D_x / max(D_x(:));
    objValueS = (clip_x-D_x);
    objValueS = -(1/(N*2))*sum(sum(objValueS.*objValueS));
    
    diff = objValueS - objValue;
    P    = exp(diff);
    R    = rand;
    
    if(P >= R)
        A = ASample;
        objValue = objValueS;

        displayCurrentState(clip, freqBest, ABest, bBest, nmean, nsig, maxV, dXImg, ampSonar)
        
        if(objValue > objValueBest)
            objValueBest = objValue;
            ABest = A;
            bBest = b;
            freqBest = freq;
        end
    end
    
     %figure(111);title(['best simulation: ' num2str(iter) ' error: ' num2str(P)]);hold on;
     D0=simulateData(freqBest,ABest,b,size(clip,1),size(clip,2),dXImg(1,:),ampSonar);
%      D0_x=D0(mask);
%      D0_x  = D0_x - mean(D0_x(:));
%      D0_x = D0_x / max(D0_x(:));
%      D0(mask)=D0_x;
     %error=(clip_x-D0_x);
      D0_x=0*clip;
      D0_x=D0.*Mask;
     error=0;
     %imagesc(D0_x);drawnow; colormap(gray);hold on;  pause(0.1);axis ij;hold on
     %xlabel(['freqBest: ' num2str(freqBest) ' ABest: ' num2str(ABest)]);hold on
   
end


if(mod(iter, dispRate) == 0)
    displayCurrentState(clip, freqBest, ABest, bBest, nmean, nsig, maxV, dXImg, ampSonar)
end
        

end


%%


function displayCurrentState(clip, freqBest, ABest, bBest, nmean, nsig, maxV, dXImg, ampSonar)
  %      figure(100);
  %      subplot(1,3,1); hold off; imagesc(clip); 
  %      D2 = computeValues(freqBest, ABest, bBest, nmean, nsig, maxV, dXImg, ampSonar);
  %      %hold on; subplot(1,2,1); plot(D2, 'r');  
  %              subplot(1,3,2); imagesc(D2);
  %      subplot(1,3,3); imagesc(D2-clip),       
end
