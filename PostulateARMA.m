function [m,Model]=PostulateARMA(ts,P)

%Variable ts is the time-series and P is the level of significance that will be used
%to for F tests associated with the model fitting.

%Please note that in the UT patch, I only have the 95% significance level
%tables and hence, whatever level you state as P, you will be running 95%
%confidence interval statistical tests (unless you use the statistical
%toolbox, in which case you will fully utilize the parameter P).

%Time-series object can be created by typing Data=iddata(ts); Then, this stucture is compatible with the
%ARMA fiting procedures in Matlab.

%Model is the model object representing the ARMA model fitted to the Data time-series object.
%It contains all necessary information about this model.
%Type set(idpoly) to see all the model properties and how you can get them.

%m is the mean value of the Data time-series calculated using the simple sample mean.

%Created by Dragan Djurdjanovic, Ann Arbor, MI.

%Last modified Feb. 26th 2009 in Austin, TX to accommodate for the lack of
%the "finv" function from the statistics toolbox of Matlab
load F_Distribution_Wrkspc %Loading the Fischer tables
ni1(length(ni1))=300;
ni2(length(ni2))=300;
%You do not need the line above if you have the statistical toolbox.



Data=iddata(ts);

%ts=get(Data,'y');
N=length(ts);
m=mean(ts);
[n1,n2]=size(ts);
ts=ts-m*ones(n1,n2);
Cycle=1;

%Initializing
CurrentModel=armax(Data,[2 1]);
n=1;
CurrentEstimationInfo=CurrentModel.EstimationInfo;
CurrentRSS=CurrentEstimationInfo.LossFcn*(N-4);


while Cycle
    n=n+1;
    OldModel=CurrentModel;
    OldRSS=CurrentRSS;
    CurrentModel=armax(Data,[2*n 2*n-1]);
    
    CurrentEstimationInfo=CurrentModel.EstimationInfo;
    CurrentRSS=CurrentEstimationInfo.LossFcn*(N-4*n);
    TestRatio=((OldRSS-CurrentRSS)/4)/(CurrentRSS/(N-4*n));
    %Control=finv(P,4,N-4*n); Had to change this!
    %If you have statistical toolbox, uncomment the line above
    %and comment the line below
    Control=FindFischer(4,N-4*n,f_95,ni1,ni2);
    [TestRatio Control;2*n 2*n-1;2*n-2 2*n-3];%this was just for debugging purposes.
    if TestRatio<Control
        Cycle=0;
        PreliminaryModel=OldModel;
        PreliminaryRSS=OldRSS;
        [TestRatio Control;2*n-1 2*n-2]; %this was just for debugging purposes.
    end
end
AR_Order=length(PreliminaryModel.a)-1;
MA_Order=length(PreliminaryModel.c)-1;


%Now check if the odd valued model is good
CurrentModel=armax(Data,[AR_Order-1 AR_Order-2]);
CurrentEstimationInfo=CurrentModel.EstimationInfo;
CurrentRSS=CurrentEstimationInfo.LossFcn*(N-2*AR_Order-4);
TestRatio=((CurrentRSS-PreliminaryRSS)/2)/(PreliminaryRSS/(N-(2*AR_Order-2)));

Control=finv(P,2,N-(2*AR_Order-2)); %Had to change this!
%If you have statistical toolbox, uncomment the line above
%and comment the line below
%Control=FindFischer(2,N-(2*AR_Order-2),f_95,ni1,ni2);

if TestRatio<Control
    PreliminaryModel=CurrentModel;
    PreliminaryRSS=CurrentRSS;
end

%Now, removing the unnecessary MA parameters.
AR_Order=length(PreliminaryModel.a)-1;
MA_Order=length(PreliminaryModel.c)-1;
CurrMA=MA_Order;
CurrentModel=PreliminaryModel;
CurrentRSS=PreliminaryRSS;

if CurrMA>1
    Cycle=1;
else
    Cycle=0;
    Model=PreliminaryModel;
    RSS=PreliminaryRSS;
end

while Cycle
    OldModel=CurrentModel;
    OldRSS=CurrentRSS;
    CurrMA=CurrMA-1;
    CurrentModel=armax(Data,[AR_Order CurrMA]);
    CurrentEstimationInfo=CurrentModel.EstimationInfo;
    CurrentRSS=CurrentEstimationInfo.LossFcn*(N-AR_Order-CurrMA-1);
    NumOfParams=AR_Order+CurrMA+1;
    TestRatio=((CurrentRSS-PreliminaryRSS)/1)/(PreliminaryRSS/(N-NumOfParams));
    
    %%%%%%%%%%%%%%%%%%
    %Control=FindFischer(1,NumOfParams,f_95,ni1,ni2); %This is a patch function I made since we do not have
    %the finv function (no statistics toolbox).
    %If you have statistical toolbox, uncomment the line below
    %and comment the line above
    Control=finv(P,1,NumOfParams);
    if TestRatio>Control
        Cycle=0;
        Model=OldModel;
        RSS=OldRSS;
    end
end %Done

%Making sure that I don't have a stupid model
if length(Model.a)<2
    Model=armax(Data,[2,1]);
end

%This is just for checking and debugging purposes - you can remove this
%part if you want to.
AR_Order=length(Model.a)-1
MA_Order=length(PreliminaryModel.c)-1