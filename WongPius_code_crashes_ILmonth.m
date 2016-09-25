%% notes
% Pius Wong
% ME383Q - Tme Series Analysis - Crash data analysis

clear all

%% import data
%get state-level data
[dataraw,txt,raw] = xlsread('dataSTATE.xls') ;  
time = dataraw(:,4);
TX_daily = dataraw(:,6);
IL_daily = dataraw(:,7);
monthindex = dataraw(:,2);

%% aggregate the data, weekly
num1 = IL_daily ;  %choose which dataset to analyze; THIS CHANGES
n1 = size(num1,1);
n = size(num1,1);
previousmonth=1;
monthsum = 0;
month = 1;
for i=1:n
    if i>1
        previousmonth = monthindex(i-1,1);
    end
    currentmonth = monthindex(i,1);
    if currentmonth==previousmonth  % if the same month
        monthsum = monthsum + num1(i,1);
    else    % if the next month
        num(month,1) = monthsum;
        month = month + 1;
        monthsum = num1(i,1);
    end
end
n = size(num,1);
time = [1:n]';% replace time function

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
plot(X(:,2),Y,'-b')
plot(regX,regY,'-k','LineWidth',2);
hold off
xlabel('Time (months)')
ylabel('Crashes')
%axis([1 n 0 2500 ])
title('Illinois monthly crashes, showing linear trend')

%% Get and graph detrended data
data = residuals;
plot(X(:,2),residuals,'b')
xlabel('Time (weeks)')
ylabel('Crashes minus linear trend')
%axis([1 n -3000 3000 ])
title('Illinois weekly crashes, detrended')

%% Autocorrelation of starting data
% estimated autocorrelation rho_hat
sum1=0; sum2=0;
for i=1:n  % i = shift variable +1 (corresponds with L in notes)
    for j=1:n-(i-1)  % j = t-index
        sum1 = data(j)*data(j+i-1) + sum1;
    end
    sum2 = data'*data;
    rho(i,1) = (1/(n-i))*sum1/( (1/n)*sum2);
end
threshold = ones(n-1,1)*2/sqrt(n);
shift = time-1;
plot(shift,rho,'r',shift,rho,'b')
xlabel('Time (weeks)')
ylabel('Estimated autocorrelation (rho)')
axis([1 2500 -10 10 ])
title('Autocorrelation function of Texas data')

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
%plot([0 n],[2/sqrt(n) 2/sqrt(n)],'-')
hold off
% RSS1 = sum(residuals(:,1).^2);  %rss of initial ARMA
% plot(residuals)


%% Graph actual vs. modeled detrended data
modeleddata = data-residuals;
plot(X(:,2),data,'b',X(:,2),modeleddata,'k')
xlabel('Time (months)')
ylabel('Detrended crash count')
%axis([1 n -50 60 ])
title('Actual vs. Modeled Illinois monthly crashes, detrended')
legend('Actual', 'Model')
strRSS = num2str(RSS1,'%1.3e\n');
str1 = ['Model RSS = ',strRSS];
text(20,15000,...
     str1); %,...
%     'FontSize',12)

spreadmetric = [Model.noisevariance; var(modeleddata); ...
    var(modeleddata)/Model.noisevariance];

%% find RSS vs. AR order
maxorder = 15; % even number
tic
RSSvsorder = zeros(maxorder/2,2);
for i=1:(maxorder/2)
    na = i*2;  % AR order
    nc = na-1;  % MA order
    ModelCHECK = armax(data,[na nc]); % sys = armax(data,[na nb nc nk])
    % phisCHECK = [ModelCHECK.a; ModelCHECK.da]';  % AR coeffs and stdevs
    % thetasCHECK = [ModelCHECK.c; ModelCHECK.dc]';
    CurrentEstimationInfo=ModelCHECK.EstimationInfo;
    RSSvsorder(i,1) = na;
    RSSvsorder(i,2) = CurrentEstimationInfo.LossFcn*(n-na-nc-1);
end
plot(RSSvsorder(:,1),RSSvsorder(:,2),'LineWidth',2)
title('Residual RSS vs. ARMA model AR order')
xlabel('AR Order')
ylabel('RSS')
toc  % time elapsed for this cell is....


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
na = 3;  % AR order
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
resid(data,ModelCHECK,'Corr',n)
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



%% Checking constant trend: get new parsimonious ARMA
for i=4:size(data,1)  %shift data to find Yt
    j = i-2;
        data2(j,1) = data(i,1) - data(i-1,1);
end
na = 3;  % AR order
nc = 3;  % MA order
Model2 = armax(data2,[na nc]); % sys = armax(data,[na nb nc nk]) 
residuals2 =  resid(data2,Model2,'Corr',n);
% RSS3 = sum(residuals3(:,1).^2)  %rss of ARMA(3,9)
phis2 = [Model2.a; Model2.da]';  % AR coeffs and stdevs
thetas2 = [Model2.c; Model2.dc]';
CurrentEstimationInfo=Model2.EstimationInfo;
RSS2 = CurrentEstimationInfo.LossFcn*(n-2*(size(phis2,1)-1)-4)


%% Checking period6 operator: get new parsimonious ARMA
for i=3:size(data,1)  %shift data to find Yt
    j = i-2;
        data4(j,1) = data(i,1) - data(i-1,1) + data(i-2,1);
end
na = 4;  % AR order
nc = 3;  % MA order
Model4 = armax(data4,[na nc]); % sys = armax(data,[na nb nc nk]) 
residuals4 =  resid(data4,Model4,'Corr',n);
% RSS3 = sum(residuals3(:,1).^2)  %rss of ARMA(3,9)
phis4 = [Model4.a; Model4.da]';  % AR coeffs and stdevs
thetas4 = [Model4.c; Model4.dc]';
CurrentEstimationInfo=Model4.EstimationInfo;
RSS4 = CurrentEstimationInfo.LossFcn*(n-2*(size(phis4,1)-1)-4)


%% visualize residuals (modify as necessary)
plot([1:n],residuals,'-',[1:n-2],residuals4)