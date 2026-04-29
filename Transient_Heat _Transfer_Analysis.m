clear; clc; close all;

%% =========================
% GIVEN DATA
% =========================
density = [7900, 100, 2700];      % kg/m^3
Cond    = [16, 0.1, 205];         % W/m-K
cp_vals = [500, 1000, 900];       % J/kg-K

T_in = 400;       % inner surface temperature (C)
T_out = 25;       % outer surface temperature (C)
T_initial = 25;   % initial temperature (C)

ri = 0.05;              % inner radius (m)
t_steel = 0.004;
t_ins   = 0.040;
t_al    = 0.003;
ro = ri + t_steel + t_ins + t_al;

Nr = 300;               % number of nodes
r = linspace(ri, ro, Nr)';
dr = r(2) - r(1);

t_final = 120;           % main simulation time (s)

r1 = ri + t_steel;
r2 = r1 + t_ins;

%% =========================
% MATERIAL PROPERTIES AT NODES
% =========================
k = zeros(Nr,1);
rho = zeros(Nr,1);
cp = zeros(Nr,1);

for i = 1:Nr
    if r(i) <= r1
        k(i) = Cond(1);
        rho(i) = density(1);
        cp(i) = cp_vals(1);
    elseif r(i) <= r2
        k(i) = Cond(2);
        rho(i) = density(2);
        cp(i) = cp_vals(2);
    else
        k(i) = Cond(3);
        rho(i) = density(3);
        cp(i) = cp_vals(3);
    end
end

alpha_node = k ./ (rho .* cp);
alpha_max = max(alpha_node);

%% =========================
% CHOOSE STABLE dt
% =========================
dt = 0.4 * dr^2 / alpha_max;
Nt = ceil(t_final / dt);

%% =========================
% MAIN TRANSIENT SOLUTION
% =========================
T = T_initial * ones(Nr,1);
Tnew = T;

T_hist = zeros(Nr,Nt);
time = zeros(1,Nt);

for n = 1:Nt
    T_hist(:,n) = T;
    time(n) = (n-1)*dt;

    for i = 2:Nr-1
        Fo = alpha_node(i) * dt / dr^2;
        Tnew(i) = T(i) + Fo * ( ...
            T(i+1) - 2*T(i) + T(i-1) + ...
            (dr/(2*r(i))) * (T(i+1) - T(i-1)) );
    end

    Tnew(1) = T_in;
    Tnew(Nr) = T_out;
    T = Tnew;
end

%% =========================
% MAIN TRANSIENT PLOT
% =========================
x = r - ri;
figure;
hold on

idx_plot = round(linspace(1, Nt, 8));
idx_plot = unique(idx_plot);

legend_entries = strings(length(idx_plot),1);

for m = 1:length(idx_plot)
    idx = idx_plot(m);
    plot(x, T_hist(:,idx), 'LineWidth', 1.5)
    legend_entries(m) = "t = " + num2str(time(idx),'%.1f') + " s";
end

xlabel('Position through wall thickness (m)')
ylabel('Temperature (C)')
title('Transient Temperature Profiles')
legend(legend_entries, 'Location','best')
grid on
hold off

%% =========================
% STABILITY COMPARISON
% =========================

% -------------------------
% Stable case
% -------------------------
dt_stable = 0.4 * dr^2 / alpha_max;
t_final_stable = 10;
Nt_stable = ceil(t_final_stable / dt_stable);

T_stable = T_initial * ones(Nr,1);
Tnew_stable = T_stable;
T_hist_stable = zeros(Nr, Nt_stable);
time_stable = zeros(1, Nt_stable);

for n = 1:Nt_stable
    T_hist_stable(:,n) = T_stable;
    time_stable(n) = (n-1)*dt_stable;

    for i = 2:Nr-1
        Fo = alpha_node(i) * dt_stable / dr^2;
        Tnew_stable(i) = T_stable(i) + Fo * ( ...
            T_stable(i+1) - 2*T_stable(i) + T_stable(i-1) + ...
            (dr/(2*r(i))) * (T_stable(i+1) - T_stable(i-1)) );
    end

    Tnew_stable(1) = T_in;
    Tnew_stable(Nr) = T_out;
    T_stable = Tnew_stable;
end

figure;
hold on

idx_plot_stable = round(linspace(1, Nt_stable, 8));
idx_plot_stable = unique(idx_plot_stable);
legend_entries_stable = strings(length(idx_plot_stable),1);

for m = 1:length(idx_plot_stable)
    idx = idx_plot_stable(m);
    plot(x, T_hist_stable(:,idx), 'LineWidth', 1.5)
    legend_entries_stable(m) = "t = " + num2str(time_stable(idx),'%.3f') + " s";
end

xlabel('Position through wall thickness (m)')
ylabel('Temperature (C)')
title('Stable Case: Fo < 0.5')
legend(legend_entries_stable, 'Location','best')
grid on
hold off

% -------------------------
% Unstable case
% -------------------------
dt_unstable = 8.0 * dr^2 / alpha_max;    % deliberately unstable
t_final_unstable = 20;                   % long enough to show instability
Nt_unstable = ceil(t_final_unstable / dt_unstable);

T_unstable = T_initial * ones(Nr,1);
Tnew_unstable = T_unstable;
T_hist_unstable = zeros(Nr, Nt_unstable);
time_unstable = zeros(1, Nt_unstable);

actual_steps_unstable = Nt_unstable;

for n = 1:Nt_unstable
    T_hist_unstable(:,n) = T_unstable;
    time_unstable(n) = (n-1)*dt_unstable;

    for i = 2:Nr-1
        Fo = alpha_node(i) * dt_unstable / dr^2;
        Tnew_unstable(i) = T_unstable(i) + Fo * ( ...
            T_unstable(i+1) - 2*T_unstable(i) + T_unstable(i-1) + ...
            (dr/(2*r(i))) * (T_unstable(i+1) - T_unstable(i-1)) );
    end

    Tnew_unstable(1) = T_in;
    Tnew_unstable(Nr) = T_out;
    T_unstable = Tnew_unstable;

    % stop if solution blows up badly
    if any(isnan(T_unstable)) || any(isinf(T_unstable)) || max(abs(T_unstable)) > 1e6
        actual_steps_unstable = n;
        break
    end
end

% trim unused columns if stopped early
T_hist_unstable = T_hist_unstable(:,1:actual_steps_unstable);
time_unstable = time_unstable(1:actual_steps_unstable);

figure;
hold on

idx_plot_unstable = round(linspace(1, actual_steps_unstable, min(8,actual_steps_unstable)));
idx_plot_unstable = unique(idx_plot_unstable);
legend_entries_unstable = strings(length(idx_plot_unstable),1);

for m = 1:length(idx_plot_unstable)
    idx = idx_plot_unstable(m);
    plot(x, T_hist_unstable(:,idx), 'LineWidth', 1.5)
    legend_entries_unstable(m) = "t = " + num2str(time_unstable(idx),'%.4f') + " s";
end

xlabel('Position through wall thickness (m)')
ylabel('Temperature (C)')
title('Unstable Case: Fo > 0.5')
legend(legend_entries_unstable, 'Location','best')
grid on
hold off

%% =========================
% SENSITIVITY TO DELTA r
% =========================
t_compare = 2;

% coarse grid
Nr_c = 20;
r_c = linspace(ri, ro, Nr_c)';
dr_c = r_c(2) - r_c(1);

k_c = zeros(Nr_c,1);
rho_c = zeros(Nr_c,1);
cp_c = zeros(Nr_c,1);

for i = 1:Nr_c
    if r_c(i) <= r1
        k_c(i) = Cond(1);
        rho_c(i) = density(1);
        cp_c(i) = cp_vals(1);
    elseif r_c(i) <= r2
        k_c(i) = Cond(2);
        rho_c(i) = density(2);
        cp_c(i) = cp_vals(2);
    else
        k_c(i) = Cond(3);
        rho_c(i) = density(3);
        cp_c(i) = cp_vals(3);
    end
end

alpha_c = k_c ./ (rho_c .* cp_c);
dt_c = 0.4 * dr_c^2 / max(alpha_c);
Nt_c = ceil(t_compare / dt_c);

T_c = T_initial * ones(Nr_c,1);
Tnew_c = T_c;

for n = 1:Nt_c
    for i = 2:Nr_c-1
        Fo = alpha_c(i) * dt_c / dr_c^2;
        Tnew_c(i) = T_c(i) + Fo * ( ...
            T_c(i+1) - 2*T_c(i) + T_c(i-1) + ...
            (dr_c/(2*r_c(i))) * (T_c(i+1) - T_c(i-1)) );
    end

    Tnew_c(1) = T_in;
    Tnew_c(Nr_c) = T_out;
    T_c = Tnew_c;
end

% fine grid
Nt_f = ceil(t_compare / dt);
T_f = T_initial * ones(Nr,1);
Tnew_f = T_f;

for n = 1:Nt_f
    for i = 2:Nr-1
        Fo = alpha_node(i) * dt / dr^2;
        Tnew_f(i) = T_f(i) + Fo * ( ...
            T_f(i+1) - 2*T_f(i) + T_f(i-1) + ...
            (dr/(2*r(i))) * (T_f(i+1) - T_f(i-1)) );
    end

    Tnew_f(1) = T_in;
    Tnew_f(Nr) = T_out;
    T_f = Tnew_f;
end

figure;
hold on
plot(r_c - ri, T_c, 'LineWidth', 1.5)
plot(r - ri, T_f, '--', 'LineWidth', 1.5)
xlabel('Position through wall thickness (m)')
ylabel('Temperature (C)')
title('Sensitivity to \Deltar at t = 2 s')
legend('Coarse grid', 'Fine grid', 'Location', 'best')
grid on
hold off

%% =========================
% SENSITIVITY TO DELTA t
% =========================
dt_small = 0.1 * dr^2 / alpha_max;
dt_large = 0.49 * dr^2 / alpha_max;

Nt_small = ceil(t_compare / dt_small);
Nt_large = ceil(t_compare / dt_large);

T_small = T_initial * ones(Nr,1);
Tnew_small = T_small;

for n = 1:Nt_small
    for i = 2:Nr-1
        Fo = alpha_node(i) * dt_small / dr^2;
        Tnew_small(i) = T_small(i) + Fo * ( ...
            T_small(i+1) - 2*T_small(i) + T_small(i-1) + ...
            (dr/(2*r(i))) * (T_small(i+1) - T_small(i-1)) );
    end

    Tnew_small(1) = T_in;
    Tnew_small(Nr) = T_out;
    T_small = Tnew_small;
end

T_large = T_initial * ones(Nr,1);
Tnew_large = T_large;

for n = 1:Nt_large
    for i = 2:Nr-1
        Fo = alpha_node(i) * dt_large / dr^2;
        Tnew_large(i) = T_large(i) + Fo * ( ...
            T_large(i+1) - 2*T_large(i) + T_large(i-1) + ...
            (dr/(2*r(i))) * (T_large(i+1) - T_large(i-1)) );
    end

    Tnew_large(1) = T_in;
    Tnew_large(Nr) = T_out;
    T_large = Tnew_large;
end

figure;
hold on
plot(x, T_small, 'LineWidth', 1.5)
plot(x, T_large, '--', 'LineWidth', 1.5)
xlabel('Position through wall thickness (m)')
ylabel('Temperature (C)')
title('Sensitivity to \Deltat at t = 2 s')
legend('Smaller \Deltat', 'Larger \Deltat', 'Location', 'best')
grid on
hold off

%% =========================
% STEADY-STATE SOLUTION
% =========================
L = 1;   % axial length

R1 = log(r1/ri) / (2*pi*Cond(1)*L);
R2 = log(r2/r1) / (2*pi*Cond(2)*L);
R3 = log(ro/r2) / (2*pi*Cond(3)*L);

Rtot = R1 + R2 + R3;
q = (T_in - T_out) / Rtot;

T_int1 = T_in  - q*R1;
T_int2 = T_int1 - q*R2;

r_ss1 = linspace(ri, r1, 200);
r_ss2 = linspace(r1, r2, 200);
r_ss3 = linspace(r2, ro, 200);

T_ss1 = T_in + (T_int1 - T_in) .* (log(r_ss1/ri) ./ log(r1/ri));
T_ss2 = T_int1 + (T_int2 - T_int1) .* (log(r_ss2/r1) ./ log(r2/r1));
T_ss3 = T_int2 + (T_out - T_int2) .* (log(r_ss3/r2) ./ log(ro/r2));

r_ss = [r_ss1, r_ss2(2:end), r_ss3(2:end)];
T_ss = [T_ss1, T_ss2(2:end), T_ss3(2:end)];

%% =========================
% VALIDATION PLOT
% =========================
figure;
hold on
plot(x, T_hist(:,end), 'LineWidth', 1.8)
plot(r_ss - ri, T_ss, '--', 'LineWidth', 1.8)
xlabel('Position through wall thickness (m)')
ylabel('Temperature (C)')
title('Validation: Numerical vs Steady-State Solution')
legend('Numerical solution at final time', 'Analytical steady-state solution', 'Location', 'best')
grid on
%% Heat Flux
T_final = T_hist(:,end);

% Local temperature gradient
dTdr = gradient(T_final, dr);

% Local radial heat flux
q_flux = -k .* dTdr;   % W/m^2

% indices for each material
idx_steel = find(r <= r1);
idx_ins   = find(r > r1 & r <= r2);
idx_al    = find(r > r2);

% average heat flux in each region
q_steel = mean(q_flux(idx_steel));
q_ins   = mean(q_flux(idx_ins));
q_al    = mean(q_flux(idx_al));

fprintf('\n===== AVERAGE HEAT FLUX IN EACH MATERIAL =====\n')
fprintf('Steel heat flux      = %.4f W/m^2\n', q_steel)
fprintf('Insulation heat flux = %.4f W/m^2\n', q_ins)
fprintf('Aluminum heat flux   = %.4f W/m^2\n', q_al)