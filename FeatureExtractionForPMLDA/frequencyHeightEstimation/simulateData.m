function [clip, dX, y] = simulateData(freq, A, b, Nrows, Ncols, x, ampSonar)

dispAll = 1;
figNum = 100;
nmean = 90;
nsig = 50;
maxV = 6.1e-3;

%Sigma s parameters, Theta_not values
dX      = x(2)-x(1);

%Sine Ripple Field Parameters
if(nargin == 0)
    freq = .5;
    A    = .3;
    b    = 0;
    
    %Image size
    Nrows   = 1000;
    Ncols   = 2000;
    
    x = [0:dX:((Ncols)*dX)-dX];  %Ground Positions
end

%Rayleigh distribution parameters
alpha = 1e2;

y = A*sin(2*pi*freq.*x+b); %+s; %Ripple Height Field
Theta_0 = atand( ampSonar./x ); %Incident Angle at every point
ThetaMax = Theta_0 + atand(2*pi*freq*A); %Determine scaling values
ThetaMin = Theta_0 - atand(2*pi*freq*A); %Determine scaling values
ThetaMaxVal = normpdf(ThetaMax, nmean, nsig);
ThetaMinVal = normpdf(ThetaMin, nmean, nsig);
D = (maxV-(1/2)*(ThetaMaxVal-ThetaMinVal)) + (ThetaMaxVal-ThetaMinVal).*(cos((2*pi*freq.*x)+b)./2);
m = mod(x- 1/(4*freq), 1/freq);
littleA = A - (m*ampSonar)./(x-m);
OcclusionMask = littleA <= y;
D = D.*OcclusionMask;
D = repmat(D, [Nrows, 1]);

%S = random('rayl',alpha,Nrows, Ncols);
%D = D.*S + randn(size(D))*.05;
%figure; imagesc(D); colormap(gray);

%m = 2*pi*freq*A*cos((2*pi*freq.*x)+b); %Slope
%ClipAngle = 90- atand( ampSonar./ (x - mod(x-1/(4*freq), 1/(freq))));

clip = D;
dX = repmat(x, [Nrows, 1]);