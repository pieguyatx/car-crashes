%% notes
% Pius Wong
% ME383Q - Tme Series Analysis - Crash data analysis

clear all

%% import data
%get state-level data
[dataraw,txt,raw] = xlsread('dataCITY.xls') ;  
time = dataraw(:,4);
TX_daily = dataraw(:,6);


%% aggregate the data, weekly
num1 = TX_daily ;  %choose which dataset to analyze; THIS CHANGES
n1 = size(num1,1);
num1 = num1(5:n1-4,1);  % start later to match gas data
n = size(num1,1);
week = 1;
weeksum = 0;
for i=1:n
    if mod(i,7)==0
        num(week,1) = weeksum;
        week = week + 1;
        weeksum = 0;
    else
        weeksum = weeksum + num1(i,1);
    end
end
n = size(num,1);
time = [1:n]';% replace time function

numFULL = num;
timeFULL = time;
nFULL = n;

%% First order polynomial fit in time
% assume y = B0 + B1xt + et

X = [ones(n,1), time] ;  % t terms
Y = num(1:n,1) ;  % y terms
betahat = ((X'*X)^-1)*(X'*Y) ;

residuals = zeros(n,1);
regY = zeros(n,1);
for i=1:n
    % find residuals
    residuals(i,1) = Y(i) - betahat(1)*X(i,1) - betahat(2)*X(i,2);
    regY(i,1) = betahat(1)*X(i,1) + betahat(2)*X(i,2);
end
SSE = residuals'*residuals;
sigma2 = (1/(n-2))*SSE; % variance of residuals


%% Graph data and trend
regXmin = min(X(:,2));
regXmax = max(X(:,2));
regX = regXmin:1:regXmax;
hold on
plot(X(:,2),Y,'-g')
plot(regX,regY,'-k','LineWidth',2);
hold off
xlabel('Time (weeks)')
ylabel('Crashes')
axis([1 n 50 300 ])
title('Austin weekly crashes (2007-2010), showing linear trend')

%% Get and graph detrended data
data = residuals;
plot(X(:,2),residuals,'g')
xlabel('Time (weeks)')
ylabel('Crashes minus linear trend')
axis([1 n -120 80 ])
title('Austin weekly crashes, detrended')


%% fit ARMA(n,m) using Dragan's code
ts = data;
P = 0.05;
[m,Model]=PostulateARMA(ts,P)
present(Model)  % use this command to see the confidence intervals
noiseVar1 = Model.noisevariance;
phis = [Model.a; Model.da]';  % AR coeffs and stdevs
thetas = [Model.c; Model.dc]';
CurrentEstimationInfo=Model.EstimationInfo;
RSS1 = CurrentEstimationInfo.LossFcn*(n-2*(size(phis,1)-1)-4)

%% analyze residuals of model
residuals = resid(data,Model,'Corr',n);  
% if "residuals" is white, autocorrs should mostly lie within CI
% given CI is 99%CI, so 1% of the autocorr values can lie out
axis([0 100 -0.5 1])

% hold on
% plot([0 n],[2/sqrt(n) 2/sqrt(n)],'-')
% hold off
% RSS1 = sum(residuals(:,1).^2);  %rss of initial ARMA
% plot(residuals)


%% Graph actual vs. modeled detrended data
modeleddata = data-residuals;
plot(X(:,2),data,'g',X(:,2),modeleddata,'k')
xlabel('Time (weeks)')
ylabel('Detrended crash count')
axis([0 n -120 80 ])
title('Actual vs. Modeled Austin weekly crashes, detrended')
legend('Actual', 'Model')
strRSS = num2str(RSS1,'%1.3e\n');
str1 = ['Model RSS = ',strRSS];
text(40,-100,...
     str1); %,...
%     'FontSize',12)

spreadmetric = [Model.noisevariance; var(modeleddata); ...
    var(modeleddata)/Model.noisevariance];



%% find & plot roots of default ARMA
polycoefficients = phis(:,1)';
r = roots(polycoefficients);
nroots = size(phis,1)-1;
rcoords = zeros(nroots,2);
rmagnitudes = zeros(nroots,1);
rangles = zeros(nroots,1);
rperiods = zeros(nroots,1);
for i=1:nroots
    magnitude = abs(r(i,1));
    angle = atan2(imag(r(i,1)),real(r(i,1)));
    rcoords(i,1) = magnitude*cos(angle); % xcoord 
    rcoords(i,2) = magnitude*sin(angle); % ycoord
    rmagnitudes(i,1) = magnitude; % root magnitude
    rangles(i,1) = angle; % root angle
    if rangles(i,1)<0
        rangles(i,1) = rangles(i,1)+2*pi();
    end
end
plot(rcoords(:,1),rcoords(:,2),'o')
axis([-1.1 1.1 -1.1 1.1])
axis square; title('Roots of AR part')
% RSS/(N-r) = Model.noisevariance;


%% fit arbitrary ARMA(n,m) (change n,m as needed to check)
na = 16;  % AR order
nc = na-1;  % MA order
ModelCHECK = armax(data,[na nc]); % sys = armax(data,[na nb nc nk])
% present(Model)  % use this command to see the confidence intervals
noiseVarCHECK = ModelCHECK.noisevariance;
% residualsCHECK = resid(data,ModelCHECK,'Corr',n);
% RSSCHECK = sum(residualsCHECK(:,1).^2)  %rss of ARMA(arbitrary)
axis([0 n -0.5 1])
hold on
plot([0 n],[2/sqrt(n) 2/sqrt(n)],'-')
hold off
phisCHECK = [ModelCHECK.a; ModelCHECK.da]';  % AR coeffs and stdevs
thetasCHECK = [ModelCHECK.c; ModelCHECK.dc]';
CurrentEstimationInfo=ModelCHECK.EstimationInfo;
RSSCHECK = CurrentEstimationInfo.LossFcn*(n-na-nc-1)
residuals = resid(data,ModelCHECK,'Corr',n)
axis([0 n -.5 1])

%% find & plot roots of selected ARMA (use cell above)
polycoefficients = phisCHECK(:,1)';
r = roots(polycoefficients);
nroots = size(phisCHECK,1)-1;
rcoords = zeros(nroots,2);
rmagnitudes = zeros(nroots,1);
rangles = zeros(nroots,1);
for i=1:nroots
    magnitude = abs(r(i,1));
    angle = atan2(imag(r(i,1)),real(r(i,1)));
    rcoords(i,1) = magnitude*cos(angle); % xcoord 
    rcoords(i,2) = magnitude*sin(angle); % ycoord
    rmagnitudes(i,1) = magnitude; % root magnitude
    rangles(i,1) = angle; % root angle
    if rangles(i,1)<0
        rangles(i,1) = rangles(i,1)+2*pi();
    end
end
plot(rcoords(:,1),rcoords(:,2),'o')
axis([-1.1 1.1 -1.1 1.1])
axis square; title('Roots of AR part')
% RSS/(N-r) = Model.noisevariance;

%% 4-step ahead prediction based on arbitrary ARMA(1,0) ONLY, code cells above
GF1 = GreenFunction(ModelCHECK,n); 
ARorder = size(phisCHECK,1)-1;
MAorder = size(thetasCHECK,1)-1;
% numFULL timeFULL nFULL
% initialize vars
prediction4 = zeros(nFULL,1);
predictionupper = prediction4;
predictionlower = prediction4;
prediction3 = prediction4;
prediction2 = prediction4;
prediction1 = prediction4;
for i=1:nFULL  % actual detrended data, based on 2004-2008 data
    Yhatpredict = betahat(1) + betahat(2)*timeFULL(i);
    actualdata = numFULL - Yhatpredict;
end
modelnoise = prediction4;
modelnoise(1:size(residuals,1),1) = residuals;
modeldata = prediction4;
modeldata(1:size(residuals,1),1) = ...
    actualdata(1:size(residuals,1),1) - residuals ;

for i=(split+1):nFULL
    % get model data
    pastdata = actualdata((i-ARorder):(i-1),1);  % past actual data
    noise0 = modelnoise((i-MAorder):i,1); % past and present noise
    modeldata(i,1) = -phisCHECK(2:ARorder+1,1)'*flipud(pastdata) + ...
        thetasCHECK(1:MAorder+1,1)'*flipud(noise0); % like xt1hat
    % get current noise / noise from 2008-2010
    modelnoise(i,1) = actualdata(i,1) - modeldata(i); 
    prediction1(i,1) = -phisCHECK(2,1)*actualdata(i,1);
    prediction2(i,1) = -phisCHECK(2,1)*prediction1(i,1);
    prediction3(i,1) = -phisCHECK(2,1)*prediction2(i,1);
    prediction4(i,1) = -phisCHECK(2,1)*prediction3(i,1);    
    % get 95% CI see p 189 of book, eq 5.3.4
    predictionerr = 1.96* sqrt(ModelCHECK.noisevariance)* ...
        sqrt(GF1(1:4)*GF1(1:4)');  % prediction interval
    predictionupper(i,1) = prediction4(i,1) + predictionerr;
    predictionlower(i,1) = prediction4(i,1) - predictionerr;
    
end

%% 4-step ahead prediction based on arbitrary ARMA(n,m), code cells above
GF1 = GreenFunction(ModelCHECK,n); 
ARorder = size(phisCHECK,1)-1;
MAorder = size(thetasCHECK,1)-1;
% numFULL timeFULL nFULL
% initialize vars
prediction4 = zeros(nFULL,1);
predictionupper = prediction4;
predictionlower = prediction4;
prediction3 = prediction4;
prediction2 = prediction4;
prediction1 = prediction4;
for i=1:nFULL  % actual detrended data, based on 2004-2008 data
    Yhatpredict = betahat(1) + betahat(2)*timeFULL(i);
    actualdata = numFULL - Yhatpredict;
end
modelnoise = prediction4;
modelnoise(1:size(residuals,1),1) = residuals;
modeldata = prediction4;
modeldata(1:size(residuals,1),1) = ...
    actualdata(1:size(residuals,1),1) - residuals ;

for i=(split+1):nFULL
    % get model data
    pastdata = actualdata((i-ARorder):(i-1),1);  % past actual data
    noise0 = modelnoise((i-MAorder):i,1); % past and present noise
    modeldata(i,1) = -phisCHECK(2:ARorder+1,1)'*flipud(pastdata) + ...
        thetasCHECK(1:MAorder+1,1)'*flipud(noise0); % like xt1hat
    % get current noise / noise from 2008-2010
    modelnoise(i,1) = actualdata(i,1) - modeldata(i); 
    % get 1-step prediction xt1hat
    pastdata1 = actualdata((i-ARorder+1):i,1);  % past+pres actual data
    noise1 = modelnoise((i-MAorder+1):i,1); % past and present noise
    prediction1(i,1) = -phisCHECK(2:ARorder+1,1)'*flipud(pastdata1) + ...
        thetasCHECK(2:MAorder+1,1)'*flipud(noise1);
    % get 2-step prediction xt2hat
    pastdata2 = actualdata((i-ARorder+2):i,1);  % past+pres actual data
    noise2 = modelnoise((i-MAorder+2):i,1); % past and present noise
    prediction2(i,1) = -phisCHECK(2,1)*prediction1(i,1) + ...
        -phisCHECK(3:ARorder+1,1)'*flipud(pastdata2) + ...
        thetasCHECK(3:MAorder+1,1)'*flipud(noise2);    
    % get 3-step prediction xt3hat
    pastdata3 = actualdata((i-ARorder+3):i,1);  % past+pres actual data
    noise3 = modelnoise((i-MAorder+3):i,1); % past and present noise
    prediction3(i,1) = -phisCHECK(2,1)*prediction2(i,1) + ...
        -phisCHECK(3,1)*prediction1(i,1) + ...
        -phisCHECK(4:ARorder+1,1)'*flipud(pastdata3) + ...
        thetasCHECK(4:MAorder+1,1)'*flipud(noise3);        
    % get 4-step prediction xt4hat
    pastdata4 = actualdata((i-ARorder+4):i,1);  % past+pres actual data
    noise4 = modelnoise((i-MAorder+4):i,1); % past and present noise
    prediction4(i,1) = -phisCHECK(2,1)*prediction3(i,1) + ...
        -phisCHECK(3,1)*prediction2(i,1) + ...
        -phisCHECK(4,1)*prediction1(i,1) + ...
        -phisCHECK(5:ARorder+1,1)'*flipud(pastdata4) + ...
        thetasCHECK(5:MAorder+1,1)'*flipud(noise4);

    % get 95% CI see p 189 of book, eq 5.3.4
    predictionerr = 1.96* sqrt(ModelCHECK.noisevariance)* ...
        sqrt(GF1(1:4)*GF1(1:4)');  % prediction interval
    predictionupper(i,1) = prediction4(i,1) + predictionerr;
    predictionlower(i,1) = prediction4(i,1) - predictionerr;
    
end

    
%% plot prediction (of detrended data)
hold on
plot(timeFULL,actualdata,'-b','LineWidth',2)
plot(timeFULL(split+1:nFULL),prediction4(split+1:nFULL),'-r','LineWidth',2)
plot(timeFULL(split+1:nFULL), predictionupper(split+1:nFULL),...
    ':k',timeFULL(split+1:nFULL),predictionlower(split+1:nFULL),':k' )
hold off
startpt = (split-50);
axis([startpt nFULL -100 100])
title('ARMA(10,9), 4-step ahead prediction; Illinois, 2008-2010')
legend('Actual','Predicted','95%CI')
xlabel('Time (week)')
ylabel('Detrended crashes');

%% plot prediction (of data)
for i=1:nFULL
    actualdata2(i,1) = betahat(1,1) + timeFULL(i,1)*betahat(2,1) ...
        + actualdata(i,1) ;
    prediction42(i,1) = betahat(1,1) + timeFULL(i,1)*betahat(2,1) ...
        + prediction4(i,1);
    predictionupper2(i,1) = prediction42(i,1) + predictionerr;
    predictionlower2(i,1) = prediction42(i,1) - predictionerr;
end
hold on
plot(timeFULL,actualdata2,'-b','LineWidth',2)
plot(timeFULL(split+1:nFULL),prediction42(split+1:nFULL),'-r','LineWidth',2)
plot(timeFULL(split+1:nFULL), predictionupper2(split+1:nFULL),...
    ':k',timeFULL(split+1:nFULL),predictionlower2(split+1:nFULL),':k' )
hold off
startpt = (split-50);
%axis([startpt nFULL -8000 10000])
title('ARMA(10,9), 4-step ahead prediction; Illinois, 2008-2010')
legend('Actual','Predicted','95%CI')
xlabel('Time (week)')
ylabel('Detrended crashes');



%% debug code
plot(timeFULL,modelnoise,'-b','LineWidth',2)
plot(timeFULL,actualdata,timeFULL,modeldata)
plot(timeFULL,modeldata,timeFULL,prediction1,timeFULL,prediction2,...
    timeFULL,prediction3,timeFULL,prediction4)
axis([startpt nFULL -3000 3000])