function [Y] = gen_correlated_sources(corr_coeff,T,Q)
Cov=ones(Q)*corr_coeff+diag(ones(Q,1)*(1-corr_coeff)); % required covariance matrix
freq=randi(100,Q,1)+1; % random frequencies between 1Hz to 100Hz
phases = 2*pi*rand(Q,1); % random phases
t=10*pi/T:10*pi/T:10*pi;
Signals=sqrt(2)*cos(2*pi*freq*t+phases); % the basic signals
if corr_coeff < 1
    A=chol(Cov)'; % Cholesky Decomposition
    Y=A*Signals;
else % Coherent Sources
    Y=repmat(Signals(1,:),Q,1);
end

