load('Data_6005311.mat');

%% Question 1
    % Manual Method
x_bar1 = mean(A1);
n = 16;
s = std(A1);
t = (x_bar1 - 0.35)/(s/sqrt(n));
p_manual = 1 - tcdf(t, 15);

    % Matlab's Built-in
% [h_q1, p_builtin, ~, stats_q1] = ttest(A1, 0.35, 'Tail', 'right');

%% Question 2
    % Manual Method
y_bar = mean(U2) - mean(A2);
num = (var(A2) / n + var(U2) / n)^2;
den = ((var(A2) / n)^2 / (n - 1)) + ((var(U2) / n)^2 / (n - 1));
DoF = num / den;

s_ybar = sqrt((var(A2)/n + (var(U2)/n)));
t_crit = tinv(1-0.05/2, DoF);
lower_CI_q2 = y_bar - t_crit * s_ybar;
upper_CI_q2 = y_bar + t_crit * s_ybar;

    % Matlab's Built-in
% [h_q2, p_builtin2, ci_q2, stats_q2] = ttest2(U2, A2, 'vartype', 'unequal')

%% Question 3
    % SW Test
[~, p_U3] = swtest(U3);
[~, p_A3] = swtest(A3);

    % Outliers Count
k = 1.4826;
% Arginine
MAD_A3 = k * median(abs(A3 - median(A3)));
threshold_A3 = 3 * MAD_A3;
sumoutliers_A3 = sum(abs(A3 - median(A3)) > threshold_A3);

% Urea
MAD_U3 = k * median(abs(U3 - median(U3)));
threshold_U3 = 3 * MAD_U3;
sumoutliers_U3 = sum(abs(U3 - median(U3)) > threshold_U3);

    % Updated SW Test
% Outliers Identification     
outliers_A3 = abs(A3 - median(A3)) > threshold_A3;
outliers_U3 = abs(U3 - median(U3)) > threshold_U3;

% Updated Data
A3_updated = A3(~outliers_A3);
U3_updated = U3(~outliers_U3);

% Updated SW Test
[~, p_A3_updated] = swtest(A3_updated);
[~, p_U3_updated] = swtest(U3_updated);

[p_q3, h_q3, STATS] = ranksum(A3_updated, U3_updated, 'method','exact');

%% Question 4
    % Control Data
[~, p_c] = ttest(C1, C2);
p_creject = sum(p_c < 0.05)/350 * 100;

    % Test Data
[~, p_t] = ttest(T1, T2);
p_nreject = sum(p_t > 0.05)/350 * 100;

    % Power
y_bar = T2(:) - T1(:); 
mu0 = 0; 
s = std(y_bar);
p0 = [mu0, s];
p1 = mean(y_bar);

power = sampsizepwr('t', p0, p1, [], 8, 'Tail','both');

    % Additional wells
n_alpha_01 = sampsizepwr('t', p0, p1, 0.9, [], 'alpha', 0.01);
n_alpha_05 = sampsizepwr('t', p0, p1, 0.9, [], 'alpha', 0.05);
n_add = n_alpha_01 - n_alpha_05;

%% Question 5
    % Manual Method
% Correlation Coefficient & Test Statistic
Sxx = sum ((Conc - mean(Conc)).^2);
Syy = sum ((OD - mean(OD)).^2);
Sxy = sum ((Conc - mean(Conc)) .* (OD - mean(OD)));
r = Sxy / sqrt(Sxx * Syy);

n = length(Conc);
t_stat = r * sqrt((n - 2) / (1 - r^2));
p_value = 2 * (1 - tcdf(abs(t_stat), n - 2));

% Coefficient of Determination
R2 = r^2;

% Regression Analysis
Beta1 = Sxy / Sxx; % Slope
Beta0 = mean(OD) - Beta1 * mean(Conc); % Intercept

Y = Beta0 + Beta1 * Conc; % Estimated model
SS_R = sum ((OD - Y).^2);
DoF = n-2;
S_epsilon = sqrt(SS_R / DoF);

S_Beta1 = (S_epsilon) / sqrt(Sxx); % Std error slope
t_crit = tinv(1- 0.05/2, DoF);
lower_CI_q5 = Beta1 - t_crit * S_Beta1;
upper_CI_q5 = Beta1 + t_crit * S_Beta1;

% Prediction Interval
SE_exp_2 = S_epsilon * sqrt(1/n + (2 - mean(Conc))^2 / Sxx);
SE_pred_2 = sqrt(SE_exp_2^2 + S_epsilon^2);

y_hat = Beta0 + Beta1 * 2;
lower_PI = y_hat - t_crit * SE_pred_2;
upper_PI = y_hat + t_crit * SE_pred_2;

    % Matlab's Built-in
% Correlation Coefficient & Test Statistic
% [r, p] = corr(Conc(:), OD(:), 'Type', 'Pearson');

% Coefficient of Determination
% R2 = r^2;

% Regression Analysis
% mdl = fitlm(Conc, OD);
% slope = mdl.Coefficients.Estimate(2);
% CI = coefCI(mdl);
% slope_CI = CI(2, :);

%Prediction Interval
% [~, PI] = predict(mdl, 2, 'Prediction', 'observation');

%% Figure 1
figure;
hold on;

% Boxplot U2 and A2
boxplot([U2, A2], 'Labels', {'Urea', 'Arginine'});

plot(1, U2, 'ro');
plot(2, A2, 'bo');
plot(1, mean(U2), 'ro', 'MarkerFaceColor', 'r');
plot(2, mean(A2), 'bo', 'MarkerFaceColor', 'b');

ylabel('Optical Density (OD)');
title('Urea vs Arginine Effect on OD (95% CIs)');
grid on;

hold off;

%% Figure 2
figure;
hold on;

% Scatter plot
scatter(1:350, -log10(p_c), 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Control Data');
scatter(1:350, -log10(p_t), 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Test Data');
yline(-log10(0.05), 'k-', 'LineWidth', 1.5, 'DisplayName', ' Significance Level, p=0.05');

xlabel('Student Index');
ylabel('p-value (-log10 scale)');
title('Control vs Test p-values');
legend('Location', 'northeast');
grid on;

hold off;

%% Figure 3
SE_exp = S_epsilon * sqrt(1/n + ((Conc - mean(Conc)).^2) / Sxx);
SE_pred = sqrt(SE_exp.^2 + S_epsilon^2);

% Confidence Bands
lower_CB = Y - t_crit * SE_exp;
upper_CB = Y + t_crit * SE_exp;

% Prediction Bands
lower_PB = Y - t_crit * SE_pred;
upper_PB = Y + t_crit * SE_pred;

figure;
hold on;

% Scatter Plot & Regression Line
scatter(Conc, OD, 'bo');
plot(Conc, Y, 'r-', 'LineWidth', 2);

% Confidence Bands
plot(Conc, lower_CB, 'k--', 'LineWidth', 1.2);
plot(Conc, upper_CB, 'k--', 'LineWidth', 1.2);

% Prediction Bands
plot(Conc, lower_PB, 'g--', 'LineWidth', 1.2);
plot(Conc, upper_PB, 'g--', 'LineWidth', 1.2);

xlabel('Antibiotic Concentration (µg/ml)');
ylabel('Optical Density (OD)');
title('Regression with Confidence and Prediction Bands');
legend('Data', 'Regression Line', 'Confidence Bands', 'Prediction Bands', 'Location', 'northeast');

textbox = sprintf('r = %.2f\nR^2 = %.2f\nSlope = %.2f', r, R2, Beta1);
text(0.1, 0.1, textbox);
grid on;

hold off;


