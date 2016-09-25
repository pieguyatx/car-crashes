%% notes
% Pius Wong
% ME383Q - Tme Series Analysis - crash test code

clear all
[num,txt,raw] = xlsread('test.xls') ;  % needs excel 97-2003

%% Part A - First order polynomial fit
% assume y = B0 + B1xt + et

n = size(num,1);
time = [1:n]';
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
sigma2 = (1/(n-2))*SSE;

regXmin = min(X(:,2));
regXmax = max(X(:,2));
regX = regXmin:1:regXmax;
plot(X(:,2),Y,'-',regX,regY,'-')
xlabel('Time (days)')
ylabel('# crashes')

data = residuals;
plot(X(:,2),residuals)
title('Detrended retail sales data')
xlabel('Time (days)')
ylabel('Crashes minus linear trend')

%% fit ARMA(n,n-1) and ARMA(n+2,n+1), compare
% see http://www.mathworks.com/help/toolbox/ident/ref/armax.html
clear phis thetas x_hat noise_hat sys 
Data = iddata(residuals);
na = 2;  % AR order
nc = na-1;  % MA order
% sys = armax(data,[na nb nc nk]) 
sys = armax(Data,[na nc])
phis = [sys.a; sys.da]'  % AR coeffs and stdevs
thetas = [sys.c; sys.dc]'    % MA coeffs and stdevs
startindex = size(phis,2);  % number of phis + 1

residuals2 = resid(sys,Data);
noise_hat = residuals2.y;
NoiseVariance = sys.noisevariance
RSS = sum(noise_hat(:,1).^2)
x_hat = data-noise_hat;
meanx = mean(x_hat)

% visualize
resid(sys,Data,'Corr',size(data,1))  % autocorrelation e and u
plot(X(:,2), data, X(:,2), x_hat)
plot(X(:,2), residuals2)


%% fit ARMA(n,m) using Dragan's code

ts = data;
P = 0.05;
[m,Model]=PostulateARMA(ts,P)

%% extra
% % residuals calculated above already
% time = num(:,3);
% timex = time(order+1:n,1);
% plot(time,residuals,'-', timex,Yhat,'-') ;  % this is "a" and "ahat"
% title('AR(4) model for residuals of linear regression')
% legend('original data','AR(4) data')
% xlabel('Time (months)')
% ylabel('Sales ($Mil)')
% 
% for i = startindex:size(data,1);     % get noise at's; col1=time, col2=data
%     % see pg 42 of book
%     noise_hat(i,1) = phis(1,1:size(phis,2))*...
%             flipud([data(i-startindex+1:i)]) ...
%        + thetas(1,2:size(thetas,2))*...
%             flipud([noise_hat(i-(size(thetas,2)-1):(i-1),1)]);
%     x_hat(i,1) = data(i) + noise_hat(i,1);
% end