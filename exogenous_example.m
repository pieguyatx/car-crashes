%% Example of VARX estimation a with VARMAX(1,1,2) process with 2 independent trends

fprintf('Example of VARX estimation a with VARMAX(1,1,2) process with 2 independent trends ...\n');

% This model is a VARMAX(1,1,2) with 2 regression parameters such that the reduced-form model is
%	Y(t) = AR*Y(t-1) + MA*W(t-1) + X(t)*b + W(t)

T = 200;	% number of samples for process
TF = 50;	% forecast period

AR = [ 1 -0.1; 0.6 0.9 ];					% initial AR(1) matrix (2 x 2 matrix)
MA = [ 0.1 0; 0.02 0.8 ];					% initial MA(1) matrix (2 x 2 matrix)
b = [ 0.001; -0.002 ];						% initial b vector (2 regression parameters)
Q = [ 0.001 0.000001; 0.000001 0.001 ];		% initial Q matrix (innovations covariance)

SpecTrue = vgxset('constant', false, 'AR', AR, 'MA', MA, 'b', b, 'Q', Q);
SpecTrue = vgxset(SpecTrue, 'Info', 'Bivariate VARMAX(1,1,2) Model');
SpecTrue = vgxset(SpecTrue, 'Series', {'Series1', 'Series2' });

% set up exogenous inputs with trends

% the model has has independent time trends
%	Y1(t) = ... + t*b1 + ...
%	Y2(t) = ... + t*b2 + ...
% which becomes
%	Y(t) = ... + t*eye(2)*b + ...

X = cell(T,1);											% exogenous inputs
XF = cell(TF,1);
for t = 1:T												% for estimation period
	X{t} = t*eye(2);
end
for t = 1:TF											% for forecast period
	XF{t} = (t + T)*eye(2);
end

Y0 = [ 1 1 ];											% initial Y(0)
n = 2;
nMA = 2;
W0 = (chol(SpecTrue.Q)'*randn(n,nMA))';					% initial W(0)

[Y, W] = vgxsim(SpecTrue, T, X, Y0, W0);				% simulate processes

% get pure VARX estimates

ESpec = vgxset('n', 2, 'constant', false, 'nAR', 1, 'nX', 2);	% set up generic VARX(1,2) model
ESpec = vgxset(ESpec, 'Info', 'Bivariate VARX(1,2) Model');
ESpec = vgxset(ESpec, 'Series', {'Series1', 'Series2' });

ESpec = vgxvarx(ESpec, Y, X, Y0);								% estimate model

% display results side-by-side

vgxdisp(SpecTrue, ESpec);

% do forecast with forecast exogenous data

[YF, YCovF] = vgxpred(ESpec, TF, XF, Y);

vgxplot(ESpec, Y, YF, YCovF);
