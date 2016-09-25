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

%% First order polynomial fit
% assume y = B0 + B1xt + et
num = TX_daily ;  %choose which dataset to analyze; THIS CHANGES

n = size(num,1);
X = [ones(n-1,1), num(2:n,1)] ;  % t terms
Y = num(2:n,1) ;  % y terms
betahat = ((X'*X)^-1)*(X'*Y) ;

residuals = zeros(n-1,1);
regY = zeros(n-1,1);
for i=1:n
    % find residuals
    residuals(i,1) = Y(i) - betahat(1)*X(i,1) - betahat(2)*X(i,2);
    regY(i,1) = betahat(1)*X(i,1) + betahat(2)*X(i,2);
end
SSE = residuals'*residuals;
sigma2 = (1/(n-1-2))*SSE;

regXmin = min(X(:,2));
regXmax = max(X(:,2));
regX = regXmin:1:regXmax;
plot(X(:,2),Y,'-',regX,regY,'-')
xlabel('Time (days)')
ylabel('Crashes')

data = residuals;
plot(X(:,2),residuals)
title('Detrended crash data')
xlabel('Time (days)')
ylabel('Crashes minus linear trend ($mil)')


%% fit ARMA(n,m) using Dragan's code

ts = data;
P = 0.05;
[m,Model]=PostulateARMA(ts,P) 
residuals = resid(data,Model,'Corr',194);
% present(Model)

%% find & plot roots
phis = [Model.a; Model.da]';  % AR coeffs and stdevs
thetas = [Model.c; Model.dc]';
polycoefficients = phis(:,1)';
r = roots(polycoefficients);
rcoords = zeros(13,2);
for i=1:13
    magnitude = abs(r(i,1));
    angle = atan2(imag(r(i,1)),real(r(i,1)));
    rcoords(i,1) = magnitude*cos(angle); % xcoord 
    rcoords(i,2) = magnitude*sin(angle); % ycoord
end
plot(rcoords(:,1),rcoords(:,2),'o')
axis([-1.5 1.5 -1.5 1.5])
axis square; title('Roots of AR part')
% RSS/(N-r) = Model.noisevariance;
RSS = sum(residuals(:,1).^2);  %rss of ARMA(13,12)

%% Checking linear trend: get new parsimonious ARMA
for i=3:size(data,1)  %shift data to find Yt
    j = i-2;
        data2(j,1) = data(i,1) - data(i-2,1);
end
na = 11;  % AR order
nc = 12;  % MA order
Model2 = armax(data2,[na nc]); % sys = armax(data,[na nb nc nk])
residuals2 =  resid(data2,Model2,'Corr',192);
RSS2 = sum(residuals2(:,1).^2)  %rss of ARMA(11,12)

%% Checking (1-b12) operator: get new parsimonious ARMA
for i=13:size(data,1)  %shift data to find Yt
    j = i-12;
        data3(j,1) = data(i,1) - data(i-12,1);
end
na = 1;  % AR order
nc = 12;  % MA order
Model3 = armax(data3,[na nc]); % sys = armax(data,[na nb nc nk]) 
residuals3 =  resid(data3,Model3,'Corr',194-12);
RSS3 = sum(residuals3(:,1).^2)  %rss of ARMA(1,12)

%% visualize residuals (modify as necessary)
plot([1:194],residuals,'-',[1:191],residuals4)