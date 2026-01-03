%% NH3_efficiency_Stoltze_full.m
% Computes NH3 equilibrium yield, efficiency, and Extended Stoltze site density
clear; close all; clc;

input_xlsx = 'Group 1-new.xlsx';
sheetname  = 'data';
output_xlsx = 'Group1_efficiency_results_Stoltze.xlsx';
R = 8.31446261815324; % J/mol·K

%% --- Read Input Data ---
if ~isfile(input_xlsx)
    error('File "%s" not found. Place it in current folder.', input_xlsx);
end
[~, sheets] = xlsfinfo(input_xlsx);
if ~any(strcmpi(sheetname, sheets))
    error('Sheet "%s" not found in %s.', sheetname, input_xlsx);
end
opts = detectImportOptions(input_xlsx, 'Sheet', sheetname);
T = readtable(input_xlsx, opts, 'Sheet', sheetname);
T.Properties.VariableNames = lower(strrep(T.Properties.VariableNames, ' ', '_'));

col_T   = find(contains(T.Properties.VariableNames, 't_k'));
col_P   = find(contains(T.Properties.VariableNames, 'p_bar'));
col_H2  = find(contains(T.Properties.VariableNames, 'h2'));
col_N2  = find(contains(T.Properties.VariableNames, 'n2'));
col_NH3 = find(contains(T.Properties.VariableNames, 'nh3'));

cases = table();
cases.T_K = T{:, col_T};
cases.P_bar = T{:, col_P};
cases.H2_frac = T{:, col_H2};
cases.N2_frac = T{:, col_N2};
cases.y_NH3_out = T{:, col_NH3};

% Normalize feed ratios
for i=1:height(cases)
    s = cases.H2_frac(i) + cases.N2_frac(i);
    if s > 1.001
        cases.H2_frac(i) = cases.H2_frac(i)/s;
        cases.N2_frac(i) = cases.N2_frac(i)/s;
    end
end

%% --- Shomate Coefficients (NIST) ---
shomate.H2  = [33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0];
shomate.N2  = [28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0];
shomate.NH3 = [19.99563, 49.77119, -15.37599, 1.921168, 0.189174, -53.25399, 203.8591, -45.944];

shomate_eval = @(coef,T) deal( ...
    coef(1)*(T/1000) + coef(2)*(T/1000)^2/2 + coef(3)*(T/1000)^3/3 + coef(4)*(T/1000)^4/4 - coef(5)/(T/1000) + coef(6), ...
    coef(1)*log(T/1000) + coef(2)*(T/1000) + coef(3)*(T/1000)^2/2 + coef(4)*(T/1000)^3/3 - coef(5)/(2*(T/1000)^2) + coef(7) ...
    );

%% --- Equilibrium Calculations ---
n = height(cases);
y_eq = NaN(n,1); xi_vals = NaN(n,1);

for i = 1:n
    T_K = cases.T_K(i);
    P_bar = cases.P_bar(i);
    P_Pa  = P_bar * 1e5;
    yH2_0 = cases.H2_frac(i);
    yN2_0 = cases.N2_frac(i);
    yNH3_0 = 0;

    % Gibbs free energy (J/mol)
    [H_H2_kJ, S_H2]   = shomate_eval(shomate.H2,  T_K);
    [H_N2_kJ, S_N2]   = shomate_eval(shomate.N2,  T_K);
    [H_NH3_kJ, S_NH3] = shomate_eval(shomate.NH3, T_K);
    G_H2  = H_H2_kJ*1000 - T_K*S_H2;
    G_N2  = H_N2_kJ*1000 - T_K*S_N2;
    G_NH3 = H_NH3_kJ*1000 - T_K*S_NH3;
    dG0 = 2*G_NH3 - (G_N2 + 3*G_H2);
    Kp = exp(-dG0/(R*T_K)); % dimensionless
    
    % Solve for xi (extent)
    fun = @(xi) ((2*xi/(1-2*xi))^2) / (((yN2_0-xi)/(1-2*xi))*((yH2_0-3*xi)/(1-2*xi))^3) - Kp*(P_Pa^2);
    xi0 = min([0.05*yN2_0, 0.01]);
    xi = fzero(fun,[0, min([yN2_0,yH2_0/3])*0.999]);
    denom = 1 - 2*xi;
    y_eq(i) = (yNH3_0 + 2*xi)/denom;
    xi_vals(i) = xi;
end

cases.y_NH3_eq = y_eq;
cases.xi = xi_vals;
cases.eta = cases.y_NH3_out ./ cases.y_NH3_eq * 100; % in %

%% --- Extended Stoltze Site Density μ_s ---
pH2_bar = cases.H2_frac .* cases.P_bar;
cases.mu_s = 1.91*(cases.T_K - 618) - 1.52*(pH2_bar - 46.25) + 78.61;
cases.mu_s_norm = (cases.mu_s - min(cases.mu_s)) / (max(cases.mu_s) - min(cases.mu_s));

%% --- Plots ---
% 1. Efficiency vs Temperature
figure('Name','eta_vs_T','NumberTitle','off');
plot(cases.T_K, cases.eta,'o-b','LineWidth',1.3);
xlabel('Temperature (K)'); ylabel('\eta (%)');
title('NH_3 synthesis efficiency vs Temperature');
grid on; saveas(gcf,'eta_vs_T.png');

% 2. Parity plot
figure('Name','Parity','NumberTitle','off');
plot(cases.y_NH3_eq, cases.y_NH3_out,'bo','MarkerFaceColor','b');
hold on; plot([0 max(cases.y_NH3_eq)],[0 max(cases.y_NH3_eq)],'k--');
xlabel('y_{NH3,eq}'); ylabel('y_{NH3,out}');
title('Parity: y_{out} vs y_{eq}');
grid on; saveas(gcf,'parity_yout_yeq.png');

% 3. Efficiency vs μ_s (Extended Stoltze)
figure('Name','eta_vs_mu_s','NumberTitle','off');
yyaxis left
plot(cases.T_K, cases.eta,'o-b','LineWidth',1.3); ylabel('\eta (%)');
yyaxis right
plot(cases.T_K, cases.mu_s_norm,'s--r','LineWidth',1.3); ylabel('Relative activity (\mu_s)');
xlabel('Temperature (K)');
title('Efficiency vs. Predicted Site Density (Extended Stoltze)');
legend('\eta','\mu_s (normalized)','Location','best'); grid on;
saveas(gcf,'eta_vs_mu_s.png');

%% --- Export Results ---
writetable(cases, output_xlsx, 'Sheet','efficiency_results');
fprintf('\nResults exported to %s\n', output_xlsx);
fprintf('Plots saved: eta_vs_T.png, parity_yout_yeq.png, eta_vs_mu_s.png\n');
