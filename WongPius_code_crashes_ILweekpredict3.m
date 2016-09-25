%% notes
% Pius Wong
% ME383Q - Tme Series Analysis - Crash data analysis, PREDICTION

clear all

%% import data
%get state-level data
[dataraw,txt,raw] = xlsread('dataSTATE.xls') ;  
time = dataraw(:,4);
TX_daily = dataraw(:,6);
IL_daily = dataraw(:,7);

%% aggregate the data, weekly
num1 = IL_daily ;  %choose which dataset to analyze; THIS CHANGES
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
split=261;  % index at which the data "splits" into two parts:
    % first part for modeling, 2nd part to test the prediction

num = numFULL(1:split);  % take data to 2008 for model
time = timeFULL(1:split);
n = size(num,1);

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
plot(X(:,2),Y,'-',regX,regY,'-')
xlabel('Time (weeks)')
ylabel('Crashes')
%axis([1 n 0 2500 ])
title('Texas weekly crashes, showing linear trend')

%% Get and graph detrended data
data = residuals;
plot(X(:,2),residuals,'r')
xlabel('Time (weeks)')
ylabel('Crashes minus linear trend')
axis([1 n -3000 3000 ])
title('Texas weekly crashes, detrended')


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
axis([0 n -0.5 1])
hold on
plot([0 n],[2/sqrt(n) 2/sqrt(n)],'-')
hold off
% RSS1 = sum(residuals(:,1).^2);  %rss of initial ARMA
% plot(residuals)


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
% na = 6;  % AR order
% nc = na-1;  % MA order
% ModelCHECK = armax(data,[na nc]); % sys = armax(data,[na nb nc nk])
% % present(Model)  % use this command to see the confidence intervals
% noiseVarCHECK = ModelCHECK.noisevariance;
% % residualsCHECK = resid(data,ModelCHECK,'Corr',n);
% % RSSCHECK = sum(residualsCHECK(:,1).^2)  %rss of ARMA(arbitrary)
% axis([0 n -0.5 1])
% hold on
% plot([0 n],[2/sqrt(n) 2/sqrt(n)],'-')
% hold off
% phisCHECK = [ModelCHECK.a; ModelCHECK.da]';  % AR coeffs and stdevs
% thetasCHECK = [ModelCHECK.c; ModelCHECK.dc]';
% CurrentEstimationInfo=ModelCHECK.EstimationInfo;
% RSSCHECK = CurrentEstimationInfo.LossFcn*(n-na-nc-1)
% residuals = resid(data,ModelCHECK,'Corr',n)
% axis([0 n -.5 1])

%% find & plot roots of selected ARMA (use cell above)
% polycoefficients = phisCHECK(:,1)';
% r = roots(polycoefficients);
% nroots = size(phisCHECK,1)-1;
% rcoords = zeros(nroots,2);
% rmagnitudes = zeros(nroots,1);
% rangles = zeros(nroots,1);
% for i=1:nroots
%     magnitude = abs(r(i,1));
%     angle = atan2(imag(r(i,1)),real(r(i,1)));
%     rcoords(i,1) = magnitude*cos(angle); % xcoord 
%     rcoords(i,2) = magnitude*sin(angle); % ycoord
%     rmagnitudes(i,1) = magnitude; % root magnitude
%     rangles(i,1) = angle; % root angle
%     if rangles(i,1)<0
%         rangles(i,1) = rangles(i,1)+2*pi();
%     end
% end
% plot(rcoords(:,1),rcoords(:,2),'o')
% axis([-1.1 1.1 -1.1 1.1])
% axis square; title('Roots of AR part')
% % RSS/(N-r) = Model.noisevariance;

    %% 5-step ahead prediction based on ARMA, code cells above
    GF1 = GreenFunction(Model,n);
    ARorder = size(phis,1)-1;
    MAorder = size(thetas,1)-1;
    % numFULL timeFULL nFULL
    % initialize vars
    prediction4 = zeros(nFULL,1);
    predictionupper = prediction4;
    predictionlower = prediction4;
    prediction3 = prediction4;
    prediction2 = prediction4;
    prediction1 = prediction4;
    prediction5 = prediction4;
    actualdata = prediction4;
    for i=1:nFULL  % actual detrended data, based on 2004-2008 data
        Yhatpredict = betahat(1) + betahat(2)*timeFULL(i);
        actualdata(i,1) = numFULL(i,1) - Yhatpredict;
    end
    modelnoise = prediction4;
    modelnoise(1:size(residuals,1),1) = residuals;
    modeldata = prediction4;
    modeldata(1:size(residuals,1)-3,1) = ... %defined up to t-4
        actualdata(1:size(residuals,1)-3,1) - residuals(1:size(residuals,1)-3,1) ;
    % i.e. 5-step ahead prediction starts from t-4, to get prediction of
    % t+1 (next future step); model & data at t-4 are known, but not after
    
    for i=(split+1):(nFULL-1)
        % get model data
%        pastdata = actualdata((i-4-ARorder):(i-5),1);
        % past actual data from only the relevant chunk of time
%        noise0 = modelnoise((i-4-MAorder):(i-4),1); % past relevant noise
        % get 1-step prediction xt1hat (3 wks ago)
        pastdata1 = actualdata((i-ARorder+1-4):(i-4),1);  % past actual data, to t-4
        noise1 = modelnoise((i-MAorder+1-4):(i-4),1); % past and present noise
        if isempty(noise1)
            noise1=0;
        end
        thetacondition2 = thetas(2:MAorder+1,1);
        if isempty(thetacondition2)
            thetacondition2 = 0;
        end
        prediction1(i-3,1) = -phis(2:ARorder+1,1)'*flipud(pastdata1) + ...
            thetacondition2'*flipud(noise1);
        % get 2-step prediction xt2hat (2 wks ago)
        pastdata2 = actualdata((i-ARorder+2-4):(i-4),1);  % past actual data, to t-4
        if isempty(pastdata2)
            pastdata2=0;
        end
        noise2 = modelnoise((i-MAorder+2-4):(i-4),1); % past noise to t-4
        if isempty(noise2)
            noise2=0;
        end
        phicondition3 = phis(3:ARorder+1,1);
        if isempty(phicondition3)
            phicondition3 = 0;
        end
        thetacondition3 = thetas(3:MAorder+1,1);
        if isempty(thetacondition3)
            thetacondition3 = 0;
        end
        prediction2(i-2,1) = -phis(2,1)*prediction1(i-3,1) + ...
            -phicondition3'*flipud(pastdata2) + ...
            thetacondition3'*flipud(noise2);
        % get 3-step prediction xt3hat (1 wk ago)
        pastdata3 = actualdata((i-ARorder+3-4):(i-4),1);  % past actual data to t-4
        if isempty(pastdata3)
            pastdata3=0;
        end
        noise3 = modelnoise((i-MAorder+3-4):(i-4),1); % past noise to t-4
        if isempty(noise3)
            noise3=0;
        end
        phicondition4 = phis(4:ARorder+1,1);
        if isempty(phicondition4)
            phicondition4 = 0;
        end
        thetacondition4 = thetas(4:MAorder+1,1);
        if isempty(thetacondition4)
            thetacondition4 = 0;
        end
        if ARorder>1
            prediction3(i-1,1) = -phis(2,1)*prediction2(i-2,1) + ...
                -phis(3,1)*prediction1(i-3,1) + ...
                -phicondition4'*flipud(pastdata3) + ...
                thetacondition4'*flipud(noise3);
        elseif ARorder==1
            prediction3(i-1,1) = -phis(2,1)*prediction2(i-2,1) + ...
                -phicondition4'*flipud(pastdata3) + ...
                thetacondition4'*flipud(noise3);
        end
        % get 4-step prediction xt4hat (current week)
        pastdata4 = actualdata((i-ARorder+4-4):(i-4),1);  % past actual data to t-4
        if isempty(pastdata4)
            pastdata4=0;
        end
        noise4 = modelnoise((i-MAorder+4-4):(i-4),1); % past noise to t-4
        if isempty(noise4)
            noise4=0;
        end
        thetacondition5 = thetas(5:MAorder+1,1);
        if isempty(thetacondition5)
            thetacondition5 = 0;
        end
        phicondition5 = phis(5:ARorder+1,1);
        if isempty(phicondition5)
            phicondition5 = 0;
        end
        if ARorder>2
            prediction4(i,1) = -phis(2,1)*prediction3(i-1,1) + ...
                -phis(3,1)*prediction2(i-2,1) + ...
                -phis(4,1)*prediction1(i-3,1) + ...
                -phicondition5'*flipud(pastdata4) + ...
                thetacondition5'*flipud(noise4);
        elseif ARorder==2
            prediction4(i,1) = -phis(2,1)*prediction3(i-1,1) + ...
                -phis(3,1)*prediction2(i-2,1) + ...
                -phicondition5'*flipud(pastdata4) + ...
                thetacondition5'*flipud(noise4);
        elseif ARorder==1
            prediction4(i,1) = -phis(2,1)*prediction3(i-1,1) + ...
                -phicondition5'*flipud(pastdata4) + ...
                thetacondition5'*flipud(noise4);
        end
        % get 5-step prediction xt5hat (next week)
        pastdata5 = actualdata((i-ARorder+5-4):(i-4),1);  % past actual data to t-4
        if isempty(pastdata5)
            pastdata5=0;
        end
        noise5 = modelnoise((i-MAorder+5-4):(i-4),1); % past noise to t-4
        if isempty(noise5)
            noise5=0;
        end        
        phicondition6 = phis(6:ARorder+1,1);
        if isempty(phicondition6)
            phicondition6 = 0;
        end
        thetacondition6 = thetas(6:MAorder+1,1);
        if isempty(thetacondition6)
            thetacondition6 = 0;
        end
        if ARorder>3
            prediction5(i+1,1) = -phis(2,1)*prediction4(i,1) + ...
                -phis(3,1)*prediction3(i-1,1) + ...
                -phis(4,1)*prediction2(i-2,1) + ...
                -phis(5,1)*prediction1(i-3,1) + ...
                -phicondition6'*flipud(pastdata5) + ...
                thetacondition6'*flipud(noise5);
        elseif ARorder==3
            prediction5(i+1,1) = -phis(2,1)*prediction4(i,1) + ...
                -phis(3,1)*prediction3(i-1,1) + ...
                -phis(4,1)*prediction2(i-2,1) + ...
                -phicondition6'*flipud(pastdata5) + ...
                thetacondition6'*flipud(noise5);
        elseif ARorder==2
            prediction5(i+1,1) = -phis(2,1)*prediction4(i,1) + ...
                -phis(3,1)*prediction3(i-1,1) + ...
                -phicondition6'*flipud(pastdata5) + ...
                thetacondition6'*flipud(noise5);
        elseif ARorder==1
            prediction5(i+1,1) = -phis(2,1)*prediction4(i,1) + ...
                -phicondition6'*flipud(pastdata5) + ...
                thetacondition6'*flipud(noise5);
        end
        % get 95% CI see p 189 of book, eq 5.3.4
        predictionerr = 1.96* sqrt(Model.noisevariance)* ...
            sqrt(GF1(1:5)*GF1(1:5)');  % for 5-step ahead
        predictionupper(i+1,1) = prediction5(i+1,1) + predictionerr;
        predictionlower(i+1,1) = prediction5(i+1,1) - predictionerr;
    end

    
%% plot prediction
hold on
plot(timeFULL,actualdata,'-b','LineWidth',2)
plot(timeFULL(split+1:nFULL),prediction5(split+1:nFULL),'-r','LineWidth',2)
plot(timeFULL(split+1:nFULL), predictionupper(split+1:nFULL),...
    ':k',timeFULL(split+1:nFULL),predictionlower(split+1:nFULL),':k' )
hold off
startpt = 1; %(split-50);
axis([startpt nFULL -5000 5000])
title('ARMA(7,6) 5-step ahead prediction; Illinois, Jun 2006-2010')
legend('Actual','Predicted','95%CI')
xlabel('Time (week)')
ylabel('Detrended crashes');

actualprederror = actualdata(split+1:size(prediction5))...
    -prediction5(split+1:size(prediction5));
varactualprederror = var(actualprederror)

%% debug code
% plot(timeFULL,modelnoise,'-b','LineWidth',2)
plot(timeFULL,actualdata,timeFULL,modeldata,timeFULL,prediction1)

plot(timeFULL,actualdata,timeFULL(split+1:size(prediction5)),...
    prediction1(split+1:size(prediction5)),...
    timeFULL(split+1:size(prediction5)),...
    prediction5(split+1:size(prediction5)))
title('1- and 5-step ahead predictions, Jun 2006-2010')
xlabel('Time (week)')
ylabel('Crashes, detrended')
legend('Actual data','1-step prediction','5-step prediction')
%axis([startpt nFULL -5000 3000])