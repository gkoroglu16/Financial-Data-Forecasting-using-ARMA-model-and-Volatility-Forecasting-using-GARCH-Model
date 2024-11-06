% tidying up
clc;    % clear the command window
clear;  % clear the workspace
close all; % close all opened figures

% Load the data
data = readtimetable('ASX200.xlsx');
head(data,5); % check
tail(data,5); % check 

% ARMA(5,5)
tp = 6; %(0,..5)
tq = 6; %(0,..5)

% creating matrix for MSE and k (nb of parameters)
MSE = zeros(tp,tq);
k = zeros(tp,tq);

% specify and estimate each ARMA(p,q) model
y = data.Adj_Returns;
T = numel(y);

for p = 0:tp-1
    for q = 0:tq-1
        model = arima(p,0,q);
        [estmd1,~,~] = estimate(model, y, 'display', 'off');
        resid = infer(estmd1,y);

        % compute nb of parameters
        k(p+1,q+1) = p+q+1;

        % MSE - mean squared errors
        MSE(p+1, q+1) = mean(resid.^2);
    end
end
MSE = reshape(MSE,tp*tq,1);
k = reshape(k,tp*tq,1);

% Construct AIC, HQIC, and BIC
aic = exp(2*k/T).*MSE; % ".*" element-wise multiplication
hqic = (log(T)).^(2*k/T).*MSE;
bic = (sqrt(T)).^(2*k/T).*MSE;

% AIC
aic = reshape (aic,tp,tq)
[r,c] = find(aic == min(min(aic)));
best_p = r-1;
best_q = c-1;
fprintf('Best model in AIC: ARMA(%d,%d)\n', best_p, best_q);

% HQIC
hqic = reshape (hqic,tp,tq)
[r,c] = find(hqic == min(min(hqic)))
best_p = r-1
best_q = c-1
fprintf('Best model in HQIC: ARMA(%d,%d)\n', best_p, best_q);

% BIC
bic = reshape (bic,tp,tq)
[r,c] = find(bic == min(min(bic)))
best_p = r-1
best_q = c-1
fprintf('Best model in BIC: ARMA(%d,%d)\n', best_p, best_q);


% GARCH(1,1) model construction
% ARMA(1,3)
tp = 2; %(0,..1)
tq = 4; %(0,..3)

for p = 0:tp-1
    for q = 0:tq-1
        model_g = arima(p,0,q);
        [estmd1,~,~] = estimate(model_g, y, 'display', 'off');
        resid_g = infer(estmd1,y);
    end
end

garch_m = garch(1,1);
fit = estimate(garch_m,resid_g);
disp(fit);

% Plot the conditional volatility

% Estimate the conditional variances (volatility squared)
variances = infer(fit, resid_g); % This gives you the conditional variances

% Calculate conditional volatility (the square root of variances)
conditionalVolatility = sqrt(variances);
disp(conditionalVolatility);    % display
mean(conditionalVolatility);    % mean

% Plot conditional volatility
figure;
plot(data.Date(1:end-1), conditionalVolatility, 'LineWidth', 1.5);
title('ASX200 Return Volatility (annualized)');
xlabel('Time');
ylabel('Volatility');
grid on;
