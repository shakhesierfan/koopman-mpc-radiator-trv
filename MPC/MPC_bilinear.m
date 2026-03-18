clear; clc;

import casadi.*
load('T0.mat')
load('solar.mat')
load('learned_matrices.mat');
lat = 51.652350; 
lon = 5.649608;


Ts_hrs = 0.25;

Np = 20;
Ndays = 2;
nt_day = 24/Ts_hrs;
Nsim = nt_day*Ndays;


T_to_K = 273.15;

Az = [62; 42; 38; 45; 93]; % Area of Each Zone
Vzone = Az.*2.5; % Volume of Each Zone

ACH = [0.4; 0.4; 0.4; 0.4; 0.4]; % Air Change Rate [1/h] 

Awall = {
    [[5.9 11.5 2.48 3.7 3.45 7.78].*2.5 Az(1) Az(1)]
    [[5.9 7.78 5.45 8.3].*2.5 Az(2) Az(2)]
    [[5.45 3.45 3.7 2.33 3.79 5.45 6].*2.5 Az(3) Az(3)]
    [[3.79 9 5.69 7.87].*2.5 Az(4) Az(4)]
    [[5.45 7.87 5.69 7.95 5.69 13.47].*2.5 Az(5) Az(5)]
};
 
orien_wall = {
        [90 180 -90 -90 0 0 NaN NaN]
        [90 180 -90 0 NaN NaN]
        [90 180 90 180 -90 -90 0 NaN NaN]
        [90 180 -90 0 NaN NaN]
        [90 180 90 -140 -90 0 NaN NaN]
};

neighbor = {
            [0 0 0 3 3 2 -1 -1] 
            [0 1 3 0 -1 -1]
            [2 1 1 0 4 5 0 -1 -1]
            [3 0 5 5 -1 -1]
            [3 4 4 0 0 0 -1 -1]
};


Aglass   =  {
            [5.9 10.3 0 0 0 0 0 0].*1 
            [4 0    0 0 0 0].*1
            [0 0    0 0 0 0 0 0 0].*1
            [0 5.3  0 0 0 0].*1
            [0 0    0 8 3.64 0 0 0].*1
};

alpwall  = {
           [12.5 12.5 12.5 7 7 7 8 8]
           [12.5 7 7 12.5 8 8]
           [7 7 7 12.5 7 7 12.5 8 8]
           [7 12.5 7 7 8 8]
           [7 7 7 12.5 12.5 12.5 8 8]
};

spc_th   = {
            [0.5 0.5 0.5 0.4 0.4 0.4 0.5 0.5]
            [0.5 0.4 0.4 0.5 0.5 0.5]
            [0.4 0.4 0.4 0.5 0.4 0.4 0.5 0.5 0.5]
            [0.4 0.5 0.4 0.4 0.5 0.5]
            [0.4 0.4 0.4 0.5 0.5 0.5 0.5 0.5]
};


rho_wall = {
            [1300 1300 1300 1900 1900 1900 2400 2400]
            [1300 1900 1900 1300 2400 2400]
            [1900 1900 1900 1300 1900 1900 1300 2400 2400]
            [1900 1300 1900 1900 2400 2400]
            [1900 1900 1900 1300 1300 1300 2400 2400]
};

cp_wall  = {
           [1000 1000 1000 1000 1000 1000 1000 1000]
           [1000 1000 1000 1000 1000 1000]
           [1000 1000 1000 1000 1000 1000 1000 1000 1000]
           [1000 1000 1000 1000 1000 1000]
           [1000 1000 1000 1000 1000 1000 1000 1000]
};

d  = {
           [0.4 0.4 0.3 0.15 0.15 0.15 0.25 0.25]
           [0.4 0.15 0.15 0.4 0.25 0.25]
           [0.15 0.15 0.15 0.15 0.15 0.15 0.4 0.25 0.25]
           [0.15 0.4  0.15 0.15 0.25 0.25]
           [0.15 0.15 0.15 0.4 0.4 0.4 0.25 0.25]
};


Uwindow = {
           [1 1 0 0 0 0 0 0]
           [1 0 0 0 0 0]
           [0 0 0 0 0 0 0 0 0]
           [0 1 0 0 0 0]
           [0 0 0 1 1 0 0 0]
};

% Heat Transfer Coefficient (Uwindow = 1/Rwindow).
gamma       = [180; 90; 0; -90; -140]; % Surface azimuth of each window [°]: 0° = south, east negative, west positive.

absorptance = {
               [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]
               [0.5 0.5 0.5 0.5 0.5 0.5]
               [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]
               [0.5 0.5 0.5 0.5 0.5 0.5]
               [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]
};

second_g    = {
               [0.1 0.1 0 0 0 0 0 0]
               [0.1 0 0 0 0 0]
               [0 0 0 0 0 0 0 0 0]
               [0 0.1 0 0 0 0]
               [0 0 0 0.1 0.1 0 0 0]
};


SHGC        = {
               [0.5 0.5 0 0 0 0 0 0]
               [0.5 0 0 0 0 0]
               [0 0 0 0 0 0 0 0 0]
               [0 0.5 0 0 0 0]
               [0 0 0 0.5 0.5 0 0 0]
};


%nzone = size(Az, 1);

n_rad_z = [1; 1; 0; 1; 4];

%n_rad = sum(n_rad_z);

Qnom = {4700
    1700
    []
    1700
    [1700 1700 1700 1700]
}; % Nominal Heat Output of Each Radiator [W]



Nhours = Ndays*24 + ceil(Np*Ts_hrs); %numel(T0_hourly);
t_hourly = 0 : Nhours;
t_fine   = 0 : Ts_hrs : t_hourly(end);
T0 = interp1(t_hourly, T0_hourly(1:Nhours+1) + T_to_K, t_fine, 'linear')';




sol = [];
idx_day = 44; % 13th Feb
for i = 1:Ndays
    tempsol = compute_sol(idx_day, lat, lon, data_st(i), gamma, Nsim, Ts_hrs);
    idx_day = idx_day + 1;
    sol = [sol tempsol];
end
sol = [sol, tempsol(:,1:Np)];


working_start = [8; 8.5; 9; 9; 9.5];       % 08:00
lunch_time    = [12.5; 12.5; 13; 13; 13];
working_end   = [17.5; 17.5; 18; 18; 18];      % 18:00

Int_Gain = compute_Int_Gain(working_start, lunch_time, working_end, Ts_hrs, Nsim, Np);



[A, Brad, BT0, Bint, Bsol, vec_idx_zone] = continuous_coupled_zones(...
 Az, Vzone, alpwall, rho_wall, spc_th, ...
 d, cp_wall, Uwindow, Awall, neighbor, Aglass, ACH, absorptance, second_g, SHGC, orien_wall, gamma);



%% ---------- SETUP

t_T0  = (0:size(T0, 1)-1) * (Ts_hrs*3600); 
t_Int = (0:size(Int_Gain, 2)-1) * (Ts_hrs*3600); 
t_sol = (0:size(sol, 2)-1) * (Ts_hrs*3600); 


%T0_fun  = @(t) interp1(t_T0,  T0,       t, 'previous', 'extrap');
%Int_fun = @(t) interp1(t_Int, Int_Gain.',t, 'previous', 'extrap').';
%Sol_fun = @(t) interp1(t_sol, sol.',     t, 'previous', 'extrap').';

% Sizes
nx    = size(A,1);
nzone = numel(vec_idx_zone);
n_rad = numel([Qnom{:}]);

n_states = nx + nzone + 2*n_rad;

Tdes    = (20.5 + T_to_K)*ones(n_rad, 1);
delta_N = (70 - 55)/log((70 - 20)/(55 - 20));

%Tdel = 25*ones(nzone, 1)./60;
Tdel_sec = 25*ones(nzone, 1)*60;


rhs_base = @(t,s,Tsup_fun,T0h,Inth,Solh) full_rhs( ...
    t, s, A, Brad, BT0, Bint, Bsol, vec_idx_zone, Qnom, Az, delta_N, ...
    Tsup_fun, T0h, Inth, Solh, Tdes, Tdel_sec );

liftP = build_lift_params(nx, nzone, Qnom, vec_idx_zone, Tdes, xmu, xsig, Az, Vzone);
liftP.delta_N = delta_N;
liftP.xmu = xmu;
liftP.xsig = xsig;
liftP.umu = umu;
liftP.usig = usig;
liftP.dmu = dmu;
liftP.dsig = dsig;





%%




% ---------- MPC Parameters ---------------

%Np = 10;

Tsupply_opt_vec = zeros(Np, Nsim);

%Tdes_opt_vec   = zeros(n_rad, Nsim);
%Tdes_opt       = zeros(n_rad, Nsim);
%dial_opt       = zeros(n_rad, Nsim);

% -----------------------------------------



H = zeros(nzone, nx);

%rows = 1:nzone;
%cols = vec_idx_zone(:).'; 

H(sub2ind(size(H), 1:nzone, vec_idx_zone(:).')) = 1;

Q = 1e-4*eye(nObs);
%QTrad = 1e-4*eye(n_rad);
%QTsen = 1e-4*eye(nzone);
%QTw   = 1e-4*eye(n_rad);

R = 1e-8*eye(nzone);

%nz = nzone;
P          = eye(nObs);
%Q_aug      = blkdiag(Q, QTrad, QTsen, QTw);
H_aug      = [zeros(nzone, 1), H, zeros(nzone, nObs - 1 - size(H,2))];




x_true(:,1) = [17*ones(nx,1) + 273.15; 20*ones(n_rad, 1) + 273.15; 17*ones(nzone,1) + 273.15; 20*ones(n_rad, 1) + 273.15];
%x_true(vec_idx_zone(1),1) = 17 + 273.15;
%x_true(vec_idx_zone(2),1) = 18+ 273.15;
%x_true(vec_idx_zone(3),1) = 14 + 273.15;
%x_true(vec_idx_zone(4),1) = 16 + 273.15;
%x_true(vec_idx_zone(5),1) = 14 + 273.15;


%x_hat = zeros(n_states, Nsim);
x_hat(:,1)  = x_true(:,1);



hour = mod((0:(Nsim+Np))'*Ts_hrs, 24);
tau = 1.5;               % smoothing (0.5h = 30 min)
smoothStep = @(x,c,tau) 0.5*(1 + tanh((x - c)/tau));

w = smoothStep(hour, working_start(1), tau) - smoothStep(hour, working_end(1), tau);
Tmin = 17 + T_to_K + (21-17)*w;
Tmax = 20 + T_to_K + (23-20)*w;







%% MPC

cost_accum = 0;

nObs = size(A_bi, 1);

mask = ones(nObs,1);
mask(badZ) = 0;
mask(1) = 1;
mask = full(mask);

time_av = zeros(Nsim,1);
err = 0;


% ---------- Lifted model sizes ----------
Hx = zeros(nzone, n_states);
Hx(sub2ind([nzone,n_states], (1:nzone)', vec_idx_zone(:))) = 1;

% Covariances in x-space
Rx = (0.2^2) * eye(nzone);
Qz = 1e-3 * eye(nObs);

Pz = 1 * eye(nObs);

for k = 1:Nsim

    t0 = tic;
    opti = casadi.Opti();

    Tsupply = opti.variable(Np,1);

    delta_1 = opti.variable(nzone-1, Np);
    delta_2 = opti.variable(nzone-1, Np);


    z_pred    = casadi.MX.zeros(size(A_bi, 1), Np+1);
    xn = (x_hat(:,k) - xmu)./xsig;
    zk = lift_poly(xn, liftP, false);
    zk(badZ) = 0;
    zk(1) = 1;
    z_pred(:,1) = (zk - zmu)./zsig;
    for j = 1:Np
        kk = k + j - 1;
        
        uk = Tsupply(j);
        dk = [T0(kk); Int_Gain(:, kk ); sol(:, kk)];

        un = (uk - umu)./usig;
        dn = (dk - dmu)./dsig;
        T0n  = dn(1); 
        Intn = dn(2:nzone+1);
        Soln = dn(nzone+2:end);
  
        zk = z_pred(:, j);
        zu_k = zk(Iu) * un;

        z_next = A_bi*zk + B_bi*un + N_bi*zu_k + BT0_bi*T0n + Bint_bi*Intn + Bsol_bi*Soln;
        z_next = z_next .* mask;
        z_next(1) = 1;
        z_pred(:,j+1) = z_next;
        
    end
    x_pred = Ck_bi*z_pred(:, 2:end);
    x_pred = x_pred.* xsig  + xmu;

    cz = 0;
    
    for ind_i = 1:numel(vec_idx_zone)
        if vec_idx_zone(ind_i) ~= vec_idx_zone(3)   % skip zone 3
            cz = cz + 1;
    
            % predicted temps for this zone across horizon (Np x 1)
            Tz = x_pred(vec_idx_zone(ind_i), :)';   % transpose to Np x 1
    
            % enforce bounds with per-zone slacks (each is 1xNp, so transpose)
            opti.subject_to( Tz >= Tmin(k+1:k+Np) - delta_1(cz,:)' );
            opti.subject_to( Tz <= Tmax(k+1:k+Np) + delta_2(cz,:)' );
        end
    end




    w_input = 1e-3;

    J_u = w_input*((Tsupply - T_to_K)'*(Tsupply - T_to_K));

    w_slack_low  = 5;
    w_slack_high = 0.5;
    d1 = delta_1(:);
    d2 = delta_2(:);
    
    J_slack = w_slack_low*(d1.'*d1) + w_slack_high*(d2.'*d2);

    J = J_slack + J_u; %+ Jreg_du; %+ J_trk;

    opti.minimize(J);

    opti.subject_to(Tsupply >= (25 + T_to_K) *ones(Np, 1));
    opti.subject_to(Tsupply <= (75 + T_to_K) *ones(Np, 1));


    opti.subject_to(delta_1(:) >= 0);
    opti.subject_to(delta_2(:) >= 0);

    if k == 1
        opti.set_initial(Tsupply, (35 + T_to_K)*ones(Np, 1));
    else
        opti.set_initial(Tsupply, [Tsupply_opt_vec(2:end, k-1);  Tsupply_opt_vec(end, k-1)] );
    end
    opti.set_initial(delta_1, zeros(nzone-1, Np));
    opti.set_initial(delta_2, zeros(nzone-1, Np));



    opts = struct;
    opts.ipopt = struct( ...
        'tol',                1e-9, ...
        'constr_viol_tol',    1e-9, ...
        'acceptable_tol',     1e-9, ...
        'max_iter',           5000, ...
        'print_level',        5 ...
    );

    opti.solver('ipopt', opts);
    solution = opti.solve();

    Tsupply_opt_vec(:, k) = solution.value(Tsupply);

    Tsupply_opt = Tsupply_opt_vec(1, k);
    time_av(k) = toc(t0);



%% Real System
    sSim = simulate_full_nonlinear([(k-1)*Ts_hrs*3600 k*Ts_hrs*3600], Tsupply_opt, x_true(:,k), rhs_base, T0(k), Int_Gain(:,k).', sol(:,k).');
    x_true(:,k+1) = sSim(:,end);
    y_meas(:,k+1) = x_true(vec_idx_zone,k+1);
    
    u = Tsupply_opt;
    y = y_meas(:,k+1);
    
    e = solution.value(x_pred(:, 1)) - x_true(:,k+1);
    err = err + sum(e(:).^2);
    
    Ju = w_input*((u - T_to_K)*(u - T_to_K));
    viol1 = max(Tmin(k+1)-y,0);
    viol2 = max(y-Tmax(k+1),0);
    viol1(3) = 0;
    viol2(3) = 0;
    
    Jviol = w_slack_low*(viol1.'*viol1) + w_slack_high*(viol2.'*viol2);
    cost_accum = cost_accum + Ju + Jviol; %+ Jreg; %+ J_trk_real;

    Topt_prev = Tsupply_opt;
    
    %% EKF
    xk_hat_phys = x_hat(:,k);
    
    dk = [T0(k); Int_Gain(:,k); sol(:,k)];
    un = (Tsupply_opt - umu)./usig;
    dn = (dk - dmu)./dsig;
    
    T0n  = dn(1);
    Intn = dn(2:nzone+1);
    Soln = dn(nzone+2:end);
    
    xk_hat_norm = (xk_hat_phys - xmu)./xsig;
    
    zk_phys = lift_poly(xk_hat_norm, liftP, false);
    zk_phys(badZ) = 0;
    zk_phys(1)    = 1;
    
    z_hat = (zk_phys - zmu)./zsig;
    z_hat = z_hat .* mask;
    z_hat(1) = 1;
    
    zu_k  = z_hat(Iu) * un;
    
    z_pred = A_bi*z_hat ...
           + B_bi*un ...
           + N_bi*zu_k ...
           + BT0_bi*T0n ...
           + Bint_bi*Intn ...
           + Bsol_bi*Soln;
    
    z_pred = z_pred .* mask;
    z_pred(1) = 1;
    x_pred_norm = Ck_bi * z_pred;
    x_pred_phys = x_pred_norm .* xsig + xmu;

    m  = numel(Iu);
    
    Ssel = zeros(m, nObs);
    Ssel(sub2ind([m,nObs], (1:m).', Iu(:))) = 1;
    
    Fk = A_bi + N_bi * (un.' * Ssel);   % (nz x nz)
    
    Pz_pred = Fk * Pz * Fk.' + Qz;
    Pz_pred = 0.5*(Pz_pred + Pz_pred.');
    
    M = diag(mask);
    Pz_pred = M*Pz_pred*M;
    Pz_pred(1,:) = 0;  Pz_pred(:,1) = 0;   % since z(1) forced to 1
    
    % --- Measurement model in z-space
    Cphys = diag(xsig) * Ck_bi;            % z -> x_phys (without +xmu)
    Hz    = Hx * Cphys;                    % z -> y
    
    y_pred_kp1 = Hx * (Cphys*z_pred + xmu);
    innov      = y_meas(:,k+1) - y_pred_kp1;
    
    S  = Hz*Pz_pred*Hz.' + Rx;
    S  = 0.5*(S + S.') + 1e-12*eye(size(S));
    
    Kz = (Pz_pred*Hz.')/S;
    
    % --- Update z mean and covariance
    z_upd = z_pred + Kz*innov;
    
    Iz = eye(size(Pz_pred));
    Pz = (Iz - Kz*Hz)*Pz_pred*(Iz - Kz*Hz).' + Kz*Rx*Kz.';
    
    Pz = 0.5*(Pz + Pz.');
    
    z_upd = z_upd .* mask;
    z_upd(1) = 1;
    
    x_upd_norm = Ck_bi*z_upd;
    x_upd_phys = x_upd_norm .* xsig + xmu;
    
    x_hat(:,k+1) = x_upd_phys;   


end



t = (0:Nsim) * Ts_hrs;


figure; hold on;
ylim([15 23.05])

C = [
    0.0000 0.4470 0.7410   % blue
    0.8500 0.3250 0.0980   % vermillion
    0.9290 0.6940 0.1250   % yellow
    0.4940 0.1840 0.5560   % purple
    0.4660 0.6740 0.1880   % green
    ];

LS = {'-','--',':','-.','-'};

% ---- Zone temperatures
plot(t, x_hat(vec_idx_zone(1),:) - T_to_K, ...
    'Color',C(1,:), 'LineStyle',LS{1}, 'LineWidth',3);
plot(t, x_hat(vec_idx_zone(2),:) - T_to_K, ...
    'Color',C(2,:), 'LineStyle',LS{2}, 'LineWidth',3);
plot(t, x_hat(vec_idx_zone(3),:) - T_to_K, ...
    'Color',C(3,:), 'LineStyle',LS{3}, 'LineWidth',3);
plot(t, x_hat(vec_idx_zone(4),:) - T_to_K, ...
    'Color',C(4,:), 'LineStyle',LS{4}, 'LineWidth',3);
plot(t, x_hat(vec_idx_zone(5),:) - T_to_K, ...
    'Color',C(5,:), 'LineStyle',LS{5}, 'LineWidth',3);

% ---- Comfort bounds
plot(t(1:end), Tmin(1:Nsim+1) - T_to_K, ...
    'k:', 'LineWidth', 1);
plot(t(1:end), Tmax(1:Nsim+1) - T_to_K, ...
    'k--', 'LineWidth', 1);

ax = gca;
xt = 0:4:t(end);
ax.XTick = xt;

t_tick_clock = datetime(0,0,0) + hours(mod(xt, 24));
ax.XTickLabel = cellstr(datestr(t_tick_clock,'HH:MM'));

xline(24,'k:','LineWidth',1.2);

xlabel('$\mathrm{Clock\ Time}$\,[{\fontsize{7}{7}\selectfont HH:MM}]', ...
       'Interpreter','latex','FontSize',12);

ylabel('$T_{\mathrm{z}}\,[^{\circ}\mathrm{C}]$', 'Interpreter','latex', 'FontSize', 14);

legend({'$T_{\mathrm{z},1}$','$T_{\mathrm{z},2}$','$T_{\mathrm{z},3}$','$T_{\mathrm{z},4}$','$T_{\mathrm{z},5}$', '$T_{\mathrm{min}}$','$T_{\mathrm{max}}$'}, ...
       'Interpreter','latex','Location','north', 'FontSize', 10.5);


grid on

figure('Color','w');%,'Units','pixels','Position',[100 100 700 520]);

left   = 0.10;
bottom = 0.14;
height = 0.78;

mainW  = 0.62;
gap    = 0.11;
solW   = 0.001;

% Main axes (temps with yyaxis)
ax1 = axes('Position',[left bottom mainW height]);
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

% ---- X ticks as clock time ----
xlim(ax1,[t(1) t(end)]);
xt = 0:4:t(end);
ax1.XTick = xt;
t_tick_clock = datetime(0,0,0) + hours(mod(xt,24));
ax1.XTickLabel = cellstr(datestr(t_tick_clock,'HH:MM'));
ax1.XTickLabelRotation = 45;

% optional vertical line at 24h
xline(ax1,24,'k:','LineWidth',1.2);

% ---- Left axis: Tsupply ----
yyaxis(ax1,'left');
h1 = plot(ax1, t(1:end-1), Tsupply_opt_vec(1,:) - T_to_K, 'LineWidth',1.5);
h1.Color = [0 0.4470 0.7410];
ax1.YColor = h1.Color;
ylim(ax1,[18 96]);
ylabel(ax1,'$T_{\mathrm{sup}}\,[^{\circ}\mathrm{C}]$','Interpreter','latex','FontSize',14);

% ---- Right axis: To ----
yyaxis(ax1,'right');
h2 = plot(ax1, t, T0(1:numel(t)) - T_to_K, '--', 'LineWidth',1.8);
h2.Color = [0.8500 0.3250 0.0980];
ax1.YColor = h2.Color;
ylabel(ax1,'$T_{\mathrm{o}}\,[^{\circ}\mathrm{C}]$','Interpreter','latex','FontSize',14);

axSol = axes('Position',ax1.Position,'Color','none');
hold(axSol,'on'); box(axSol,'off');

axSol.XLim  = ax1.XLim;
axSol.XTick = [];          % no duplicate x-axis
axSol.XColor = 'none';

h3 = plot(axSol, t, sol(1,1:numel(t)), ':', 'LineWidth',1.6);
h3.Color = [0.4660 0.6740 0.1880];

yl = ylim(axSol);
axSol.YLim   = yl;
axSol.YTick  = [];
axSol.YColor = 'none';

xR = left + mainW + gap;
xR = min(xR, 0.98-solW); 

axSolR = axes('Position',[xR bottom solW height], ...
              'Color','none','YAxisLocation','right','Box','off');
axSolR.XTick  = [];
axSolR.XColor = 'none';
axSolR.YLim   = yl;
axSolR.YColor = h3.Color;

ylabel(axSolR,'$q_{\mathrm{sol}}\ (\mathrm{north\ walls})\,[\mathrm{W/m^2}]$', ...
    'Interpreter','latex','FontSize',14);

linkaxes([ax1 axSol],'x');

legend(ax1,[h1 h2 h3], ...
    '$T_{\mathrm{sup}}$','$T_{\mathrm{o}}$','$q_{\mathrm{sol}}$', ...
    'Interpreter','latex','Location','north','FontSize',13);

set(gcf,'PaperPositionMode','auto');
%exportgraphics(gcf,'myplot.png','Resolution',300);









function dy = ekf_bi_aug(k, A_bi, Brad_bi, BT0_bi, Bint_bi, Bsol_bi, ...
    vec_idx_zone, Qnom, Az, delta_N, ...
    Tsup_fun, T0_fun, Int_fun, Sol_fun, ...
    Tdes, Tdel_sec, ...
    Qc)

    
end

function F = jacobian_fd(fun, x)
    n = numel(x);
    fx = fun(x);
    F = zeros(n,n);

    % scale-aware step
    eps0 = 1e-6;
    for i = 1:n
        h = eps0*(1 + abs(x(i)));
        xp = x; xm = x;
        xp(i) = xp(i) + h;
        xm(i) = xm(i) - h;
        F(:,i) = (fun(xp) - fun(xm)) / (2*h);
    end
end



%function H = H_select(idx, n)
%m = numel(idx);
%H = zeros(m,n);
%for i = 1:m
%    H(i, idx(i)) = 1;
%end
%end



%% Computing External Disturbances

function sSim = simulate_full_nonlinear(tgrid, u, s0, rhs_full, T0, Int, Sol)
% Fully ZOH simulation:
%   u, T0, IntG, Sol are held constant on each [kTs,(k+1)Ts] interval.

    opts = odeset('RelTol',1e-7,'AbsTol',1e-12);

    sSim = s0;
    Ts = tgrid(2) - tgrid(1);

    for k = 1:(numel(tgrid)-1)

        % --- Freeze inputs at the LEFT endpoint (ZOH) ---
        tk   = tgrid(k);
        uk   = u(k);

        T0k   = T0(k);
        Intk  = Int(k,:).';
        Solk  = Sol(k,:).';

        % Constant functions over this interval
        Tsup_fun_k = @(t) uk;

        % Build an RHS that uses constants (no time-varying disturbance calls)
        rhs_k = @(t,s) full_rhs_zoh( ...
            s, Tsup_fun_k, T0k, Intk, Solk, rhs_full);

        % Integrate one step
        [~, sTemp] = ode45(rhs_k, [0 Ts], s0, opts);  % local time 0..Ts
        s0 = sTemp(end,:).';
        sSim = [sSim s0];
    end
    
    if ~isreal(sSim)
        error('simulate_full_nonlinear: complex states detected');
    end
end

function ds = full_rhs_zoh(s, Tsup_fun, T0k, Intk, Solk, rhs_full)

    T0_fun_k  = @(t) T0k;
    Int_fun_k = @(t) Intk;
    Sol_fun_k = @(t) Solk;

    ds = rhs_full(0, s, Tsup_fun, T0_fun_k, Int_fun_k, Sol_fun_k);
end


function ds = full_rhs(t,s, ...
    A, Brad, BT0, Bint, Bsol, ...
    vec_idx_zone, Qnom, Az, delta_N, ...
    Tsup_fun, T0_fun, Int_fun, Sol_fun, ...
    Tdes, Tdel_sec)

    nx    = size(A,1);
    nzone = numel(vec_idx_zone);
    n_rad = numel([Qnom{:}]);

    x    = s(1:nx);
    Trad = s(nx+1:nx+n_rad);
    Tsen = s(nx+n_rad+1:nx+n_rad+nzone);
    Tw   = s(nx+n_rad+nzone+1:end);

    Tsupply = Tsup_fun(t);
    T0      = T0_fun(t);
    IntG    = Int_fun(t);
    Sol     = Sol_fun(t);

    [dTrad, dTsen, dTw] = compute_T_rad_ct(x, vec_idx_zone, Tdes, Tsupply, Qnom, Trad, Tsen, Tdel_sec, delta_N, Tw);

    q_rad = compute_q_rad_ct(Qnom, Trad, x, Az, delta_N, vec_idx_zone, Tsupply);

    dx = A*x + BT0*T0 + Bint*IntG + Bsol*Sol + Brad*q_rad;

    ds = [dx; dTrad; dTsen; dTw];
end

function q_rad = compute_q_rad_ct(Qnom, Trad, x, Az, delta_N, vec_idx_zone, Tsupply)
    nzone = numel(Qnom);

    if isa(x, 'sym') == 1
        q_rad = sym(zeros(nzone,1));
    else
        q_rad = zeros(nzone,1);
    end

    count_rad = 0;
    %safe_lmtd_sym = @(d1,d2) lmtd_stable(d1,d2);

    for z = 1:nzone
        for r = 1:size(Qnom{z},2)
            count_rad = count_rad + 1;

            Tz = x(vec_idx_zone(z));

            %d1 = Tsupply - Tz;
            %d2 = 2*Trad(count_rad) - Tsupply - Tz;
            %delta = safe_lmtd(d1, d2);

            %if isa(x, 'sym') == 1
            %    delta = lmtd_stable(d1,d2);
            %else
            epsd = 1e-4;
            delta = 0.5*((Trad(count_rad) - Tz) + sqrt((Trad(count_rad) - Tz).^2 + epsd^2));       
            
            %safe_lmtd(d1, d2);
            %end


            q_rad(z) = q_rad(z) + Qnom{z}(r) * (delta/delta_N)^1.3 / Az(z);
        end
    end
end

function [dTrad, dTsen, dTw] = compute_T_rad_ct(x, vec_idx_zone, Tdes, Tsupply, Qnom, Trad, Tsen, Tdel_sec, delta_N, Tw)
nzone = numel(vec_idx_zone);
n_rad = numel([Qnom{:}]);

    if isa(x, 'sym')
        rho_w = sym(983); g = sym(9.81); a = sym(5.2965); b = sym(3.3160);
        cp = sym(4180); Crad = sym(1*10000); kappa = sym(5); U_in = sym(1000);
        Cw = sym(21*1000);
        A_rad = sym(3);
        dTrad = sym(zeros(n_rad,1));
        dTw   = sym(zeros(n_rad,1));
        dTsen = sym(zeros(nzone,1));
        R_hyd = sym(zeros(n_rad,1));
    else
        rho_w = 983; g = 9.81; a = 5.2965; b = 3.3160;
        cp = 4180; Crad = 1*10000; kappa = 5;
        Cw = 21*1000;
        U_in = 1000; 
        A_rad = 3;
        dTrad = zeros(n_rad,1);
        dTw   = zeros(n_rad,1);
        dTsen = zeros(nzone,1);
        R_hyd = zeros(n_rad,1);
    end
    
    for z=1:nzone
        dTsen(z) = (x(vec_idx_zone(z)) - Tsen(z)) / Tdel_sec(z);

        %if (Tsen(z) >= Tdes(1) + 2)
        %    disp("end")
        %end
    end
    
    counter = 0;
    for z=1:nzone
        for r=1:size(Qnom{z},2)
            counter = counter + 1;
    
            
            temp = (Tdes(counter) - Tsen(z) + 2)/4;
            
            %if Tdes(counter) - Tsen(z) <= -2
            %    disp('end');
            %end
            valve = 0.5*(1 + tanh(kappa*(temp - 0.5)));
            LS         = rho_w * 1.6 / (3600 * sqrt(1e5));

            k_v = rho_w * 0.73 / (3600 * sqrt(1e5)); %* valve;
            temp_kv = 1/(valve*k_v + 1e-8)^2  + 1/LS^2;
            R_hyd(counter) = 1/sqrt(temp_kv);
        end
    end
    
    R_sum = sum(R_hyd);
    counter = 0;
    
    for z=1:nzone
        for r=1:size(Qnom{z},2)
            counter = counter + 1;
    
            mm = sqrt((rho_w*g*a) / (1 + rho_w*g*b*R_sum^2)) * R_hyd(counter);
    
            Tz = x(vec_idx_zone(z));
            %d1 = Tsupply - Tz;
            %d2 = 2*Trad(counter) - Tsupply - Tz;
    
            %delta = safe_lmtd(d1, d2);
            %epsd = 0.2;
            %delta = 1/epsd * log1p( exp((Trad(counter) - Tz)*epsd) ) ;
            epsd = 1e-4;                 % smoothing width (K)
            delta = 0.5*((Trad(counter) - Tz) + sqrt((Trad(counter) - Tz).^2 + epsd^2));
    
            dTrad(counter) = ( U_in*A_rad*(Tw(counter) - Trad(counter)) ...
                             - Qnom{z}(r)*(delta/delta_N)^1.3 ) / Crad;

            dTw(counter) = ( 2*mm*cp*(Tsupply - Tw(counter)) - U_in*A_rad*(Tw(counter) - Trad(counter)) )/Cw;
        end
    end

end


%% ----------------------------------------------------------------------
function [A, Bu, BT0, Bint, Bsol, vec_idx_zone] = continuous_coupled_zones(...
    Az, Vzone, alpwall, rho_w, spc_th, d, cp_wall, Uwindow, Awall, neighbor, Aglass, ACH, absorptance, second_g, SHGC, orien_wall, gamma)
    
    %Constants.rho_air = 
    rho = 1.2041; %Constants.rho_air; %1.2041
    cp  = 1012; %Constants.C_air;   %1012

    rhoCp = rho * cp;

    n_wall = numel([Awall{:}]) - ( nnz([neighbor{:}] > 0)  /2   );

    n_zone = size(Az, 1);

    A    = zeros(n_wall + n_zone, n_wall + n_zone);
    BT0  = zeros(n_wall + n_zone, 1);
    Bsol = zeros(n_wall + n_zone, size(gamma, 1));
    Bu   = zeros(n_wall + n_zone, n_zone);
    Bint = zeros(n_wall + n_zone, n_zone);
    
    vec_idx_wall = cell(n_zone, 1);
    vec_idx_wall{1} = [2 3 4 5 6 7 8 9];

    vec_idx_zone = zeros(1, n_zone);
    vec_idx_zone(1) = 1;
    vec_idx_zone(2) = 10;

    for idx_zone = 2:n_zone

        
        
        %idx_temp = zeros(1, 4);
        idx_temp = zeros(1, numel(Awall{idx_zone}));

        larg_idx_wall = max([vec_idx_wall{:}]) + 1;

        neighbor_temp = neighbor;

        for i = 1:numel(Awall{idx_zone})
            if neighbor_temp{idx_zone}(i) ~= 0 && neighbor_temp{idx_zone}(i) ~= -1 && neighbor_temp{idx_zone}(i) < idx_zone 
                
                nb = neighbor_temp{idx_zone}(i);
                idx_test = find(neighbor_temp{nb}(:) == idx_zone, 1, 'last');   
                idx_temp(i) = vec_idx_wall{nb}(idx_test);

                neighbor_temp{nb}(idx_test) = -9;
                

            else
                idx_temp(i) = larg_idx_wall + 1;
            end

            larg_idx_wall = max(larg_idx_wall(:), max(idx_temp(:)));
        end

        vec_idx_wall{idx_zone} = idx_temp;
        
        if idx_zone < n_zone
            vec_idx_zone(idx_zone + 1) = max(idx_temp)+1;
        end

    end

    
    for i = 1:n_zone % Zone
        
        idx_zone = vec_idx_zone(i);

        Cz    = rhoCp * Vzone(i);
        Ginf  = rhoCp * Vzone(i) * ACH(i) / 3600;

        A(idx_zone, idx_zone) =  -Ginf/Cz;
        Bu(idx_zone, i)   = Az(i)/Cz;
        Bint(idx_zone, i) = Az(i)/Cz;
        BT0(idx_zone)     = BT0(idx_zone) + Ginf /Cz;

        

        for j = 1:numel(Awall{i})  % Wall
            
            idx_wall = vec_idx_wall{i}(j);
        
            Aopaque = Awall{i}(j) - Aglass{i}(j); %- Aframe(i, j);
            Awindow = Aglass{i}(j); %+ Aframe(i, j);

            rhoCp_w = rho_w{i}(j) * cp_wall{i}(j);
            Cwall   = rhoCp_w * Aopaque * d{i}(j);
    
            Rwall   = 1/alpwall{i}(j) + spc_th{i}(j)*d{i}(j)/2; 
            Rwindow = 1/Uwindow{i}(j);  %1/(Uwindow(i, j)*Awindow);

            % T_z

            A(idx_zone, idx_zone) = A(idx_zone, idx_zone) - Aopaque/(Rwall * Cz);
            A(idx_zone, idx_zone) = A(idx_zone, idx_zone) - Awindow/(Rwindow * Cz);
            A(idx_zone, idx_wall) = A(idx_zone, idx_wall) + Aopaque/(Rwall * Cz);

            % T_w
            
            A(idx_wall, idx_wall) = A(idx_wall, idx_wall) - Aopaque/(Rwall * Cwall);
            A(idx_wall, idx_zone) = A(idx_wall, idx_zone) + Aopaque/(Rwall * Cwall);
            
            % off diagonal
                            
            if neighbor{i}(j) == 0

                A(idx_wall, idx_wall) = A(idx_wall, idx_wall) - Aopaque/(Rwall * Cwall);

                BT0(idx_wall) = BT0(idx_wall) + Aopaque/(Rwall*Cwall);

                idx_orien = find(gamma == orien_wall{i}(j));
                if isempty(idx_orien)
                    error("Orientation %g not in gamma list.", orien_wall{i}(j));
                end
                if Awindow ~= 0
                    Bsol(idx_zone, idx_orien) = Bsol(idx_zone, idx_orien) ...
                        + second_g{i}(j)*SHGC{i}(j)*Awindow/Cz;
                end


                Bsol(idx_wall, idx_orien)  = Bsol(idx_wall, idx_orien) + absorptance{i}(j)*Aopaque/Cwall;

            elseif neighbor{i}(j) == -1

                %if j == numel(Awall{i})
                %    idx_orien = size(gamma, 1) + 1;
                %end
                %Bsol(idx_wall, idx_orien)  = Bsol(idx_wall, idx_orien) + absorptance{i}(j)*Aopaque/Cwall;

            end

            A_total = sum( Awall{i}(:) - Aglass{i}(:)); %- Aframe(i, :) );

            for idx_window = 1:size(Aglass{i}, 2)
                
                if Aglass{i}(idx_window) ~= 0
                    idx_orien = find(gamma == orien_wall{i}(idx_window));

    
                    Awindow_temp = Aglass{i}(idx_window); %+ Aframe(i, idx_window);
                    Bsol(idx_wall, idx_orien) = Bsol(idx_wall, idx_orien) ...
                        + (1 - second_g{i}(idx_window))*SHGC{i}(idx_window)*Awindow_temp*Aopaque/(Cwall*A_total);
                end
            end

                    

            % BT0
            BT0(idx_zone)  = BT0(idx_zone) + Awindow/(Cz * Rwindow);
    

            
        end

    end


    %Ad   = expm(A * Ts_hrs * 3600);

    %coef = A\(expm(A * Ts_hrs * 3600) - eye(size(A)));   

    %Bu_d    = coef*Bu; 
    %BT0_d   = coef*BT0;
    %Bint_d  = coef*Bint;
    %Bsol_d  = coef*Bsol;

end


%% Computing External Disturbances

function Int_Gain = compute_Int_Gain(working_start, lunch_time, working_end, Ts_hrs, Nsim, Np)
    
    nzone = size(lunch_time, 1);
    Ndays = ceil(Nsim/(24/Ts_hrs));

    Int_Gain_day = zeros(nzone, 24/Ts_hrs);
    
    start_idx = round(working_start / Ts_hrs) + 1;
    lunch_idx = round(lunch_time    / Ts_hrs) + 1;
    end_idx   = round(working_end   / Ts_hrs) + 1;
    
    for i = 1:nzone
        s = start_idx(i);
        l = lunch_idx(i);
        e = end_idx(i);
    
        % Morning: start -> lunch
        N1 = l - s;
        t1 = 0:N1;
        Int_Gain_day(i, s:l) = round(10 * sin((t1/N1) * pi));
    
        % Afternoon: lunch -> end
        N2 = e - l;
        t2 = 0:N2;
        Int_Gain_day(i, l:e) = round(10 * sin((t2/N2) * pi));
    end
    
    Int_Gain = repmat(Int_Gain_day, 1, Ndays);
    Int_Gain = [Int_Gain, Int_Gain_day(:,1:Np)];
end


function sol_day = compute_sol(n, lat, lon, data_sol, gamma, Nsim, Ts_hrs)
    
    Ndays = ceil(Nsim/(24/Ts_hrs));

    tz  = 1;
    %n   = 306;  % the day of the year
    beta   = 90; % the angle between the plane of the surface in question and the horizontal
    rho_g  = 0.50;
    
    phi  = deg2rad(lat);
    bet  = deg2rad(beta);
    
    B   = deg2rad((n-1)*360/365);
    
    %delta = deg2rad(23.45 * sind(360*(284 + n)/365));
    
    delta = (0.006918 - 0.399912*cos(B) + 0.070257*sin(B) - ...
        0.006758*cos(2*B) +  0.000907*sin(2*B) ...
        - 0.002697*cos(3*B) + 0.00148*sin(3*B));
    
    EoT = 229.2*(0.000075 + 0.001868*cos(B) - 0.032077*sin(B) ...
                - 0.014615*cos(2*B) - 0.04089*sin(2*B));
    
    L_st = 15*tz;
    TC   = 4*((360 - L_st) - (360 - lon)) + EoT;
    
    I = zeros(size(gamma, 1), 24);
    
    ghi = zeros(24, 1);
    dhi = zeros(24, 1);
    dni = zeros(24, 1);
    
    for j = 1:size(gamma, 1) %+ 1

        if j == size(gamma, 1) + 1
            gam  = 0;
            bet  = 0;
        else
            gam  = deg2rad(gamma(j));
        end
        for i = 1:24
            t_mid = (i - 0.5);             
            LST   = t_mid + TC/60;          
            omega = deg2rad(15*(LST - 12)); 
        
            % cos(theta_z)
            cos_tz = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(omega);
            cos_tz = max(cos_tz, 0);
        
            % cos(theta) on tilted plane (Duffie Eq. 1.6.x form)
            cos_th =  sin(delta)*sin(phi)*cos(bet) ...
                    - sin(delta)*cos(phi)*sin(bet)*cos(gam) ...
                    + cos(delta)*cos(phi)*cos(bet)*cos(omega) ...
                    + cos(delta)*sin(phi)*sin(bet)*cos(gam)*cos(omega) ...
                    + cos(delta)*sin(bet)*sin(gam)*sin(omega);
            cos_th = max(cos_th, 0);
        
            %Rb = 0;
            %if cos_tz > 0
            %    Rb = cos_th / max(cos_tz, eps_cz);
            %end
        
            ghi(i) = data_sol.intervals(i).irradiation.cloudy_sky.ghi;
            dhi(i) = data_sol.intervals(i).irradiation.cloudy_sky.dhi;
            dni(i) = data_sol.intervals(i).irradiation.cloudy_sky.dni;
        
            Fd = (1 + cos(bet))/2;
            Fg = (1 - cos(bet))/2;
        
            I_beam = dni(i) * cos_th;
            %I_beam = (ghi(i) - dhi(i)) * Rb;
            I_diff = dhi(i) * Fd;
            I_gnd  = rho_g  * ghi(i) * Fg;
        
            I(j, i) = I_beam + I_diff + I_gnd;
        
            %test(i) = dhi(i) + dni(i)*cos_tz;
        end
    
    end
    
    sol_day = repelem(I, 1, floor(1/Ts_hrs));
    
    %sol   = repmat(sol_day, 1, Ndays);
    %sol = [sol, sol_day(:,1:Np)];




end





 function Z = lift_poly(sK, P, linear)
    rho_w = 983; g = 9.81; a = 5.2965; b = 3.3160; % a and b are for pump
    cp = 4180; C_rad = 10000; Cw = 21*1000;
    A_rad = 3;
    U_in = 1000;
    rho = 1.2041;
    cp_air  = 1012;

    rhoCp = rho * cp_air;

    [~, N] = size(sK);

    if linear
        Z = [ones(1,N);
             sK];
        return;
    end

    % ---- split normalized blocks
    xK    = sK(P.idx_x,    :);
    TradK = sK(P.idx_Trad, :);
    TsenK = sK(P.idx_Tsen, :);
    TwK   = sK(P.idx_Tw, :);

    % ---- reconstruct physical (Kelvin)
    % Sensor temps
    mu_Tsen  = P.xmu(P.idx_Tsen);
    sig_Tsen = P.xsig(P.idx_Tsen);
    Tsen_phys = mu_Tsen + sig_Tsen .* TsenK;          % [nzone x N]

    % Radiator temps
    mu_Trad  = P.xmu(P.idx_Trad);
    sig_Trad = P.xsig(P.idx_Trad);
    Trad_phys = mu_Trad + sig_Trad .* TradK;


    mu_Tw  = P.xmu(P.idx_Tw);
    sig_Tw = P.xsig(P.idx_Tw);
    Tw_phys   = mu_Tw    + sig_Tw .* TwK;

    % Zone temps Tz from x(vec_idx_zone)
    Tz_phys = zeros(P.nzone, N);
    for z = 1:P.nzone
        idx = P.vec_idx_zone(z);
        Tz_phys(z,:) = P.xmu(idx) + P.xsig(idx) .* xK(idx,:);
    end

    % ---- per-radiator features
    nrad = P.nrad;
    e_phys    = zeros(nrad, N);   % thermostat error [K]
    valveK = zeros(nrad, N);             % valve opening in [0,1]
    dT_rz = zeros(nrad, N);
    dT_wr = zeros(nrad, N);
    dT_rw = zeros(nrad, N);
    delta_star = zeros(nrad, N);
    phi_star        = zeros(nrad, N);
    phi_star_zone   = zeros(P.nzone - 1, N);
    idx_room = 0;
    for r = 1:nrad
        z = P.rad_zone(r);
        Cz    = rhoCp * P.Vzone(z);

        e_phys(r,:) = P.Tdes_phys(r) - Tsen_phys(z,:);
        temp        = (e_phys(r,:) + 2)/4;
        valve      = 0.5*(1 + tanh(5*(temp - 0.5)));
        valve_temp = rho_w * 0.73 / (3600 * sqrt(1e5));
        LS         = rho_w * 1.6 / (3600 * sqrt(1e5));

        valve_temp2 = 1./(valve.*valve_temp + 1e-8).^2 + 1/LS.^2;  %1./sqrt(valve_temp + 1e-6) + 1/sqrt(1.6);

        valveK(r,:) = 1./sqrt(valve_temp2);

%U_in*A_rad*(Tw(counter) - Trad(counter)
        dT_wr(r,:)  = U_in.*A_rad.*(Tw_phys(r,:) -   Trad_phys(r,:))./ C_rad; %   
        dT_rw(r,:)  = U_in.*A_rad.* (Tw_phys(r,:) -   Trad_phys(r,:)) ./ Cw; % 

        dT_rz(r,:)  = Trad_phys(r,:) - Tz_phys(z,:);
        epsd = 1e-4;
        delta_star(r,:) = 0.5*( dT_rz(r,:)  + sqrt(dT_rz(r,:).^2 + epsd^2) );       

        %
        phi_star(r,:)   = P.Qnom(r).*((delta_star(r,:) ./ P.delta_N) .^ 1.3)./ C_rad;
        if z ~= 3
            if z > 3
                idx_room = z - 1;
            else
                idx_room = z;
            end
            phi_star_zone(idx_room, :) = phi_star_zone(idx_room, :) + P.Qnom(r).*((delta_star(r,:) ./ P.delta_N) .^ 1.3).*P.Az(z)./ Cz; ; 
        end


    end

    sum_valve = sum(valveK, 1);        % 1 x N
    RK = zeros(nrad, N); 
    R_u = zeros(nrad, N);
    for r = 1:nrad
        RK(r,:)  =  2*cp.* sqrt((rho_w*g*a) ./ (1 + rho_w*g*b*sum_valve.^2)).*valveK(r,:).*Tw_phys(r,:)./Cw; % 
        R_u(r,:) =  2*cp.* sqrt((rho_w*g*a) ./ (1 + rho_w*g*b*sum_valve.^2)).*valveK(r,:)./Cw; % 
    end

    Z = [ones(1,N); 
         sK];
        nTw = P.idx_Tw(end);
        
        valid = (1:nTw) > P.idx_x(end) | ismember(1:nTw, P.vec_idx_zone(:).');
        idx = find(valid);
        idx = idx(:);     % important
        m = numel(idx);
        
        [A,B] = ndgrid(1:m, 1:m);
        keep2 = A <= B;
        pairs = [idx(A(keep2)) idx(B(keep2))];
        
        [A,B,C] = ndgrid(m - nrad:m, m - nrad:m, m - nrad:m);
        keep3 = (A <= B) & (B <= C);
        triples = [idx(A(keep3)) idx(B(keep3)) idx(C(keep3))];
        Zp = sK(pairs(:,1),:) .* sK(pairs(:,2),:);
        Zt = sK(triples(:,1),:) .* sK(triples(:,2),:) .* sK(triples(:,3),:);
        
        Z = [Z; Zp; Zt];


    Z = [Z;
         dT_rz;
         dT_wr;
         dT_rw;
         valveK;
         sum_valve;
         phi_star_zone;
         phi_star;
         RK;
         R_u];
         %f_ed];
         %f_e2;
         %f_d2];
         %g1;
         %g2;
         %g3];

end




function P = build_lift_params(nx, nzone, Qnom, vec_idx_zone, Tdes, xmu, xsig, Az, Vzone)
    P.Az = Az;
    P.Vzone = Vzone;
    % radiator -> zone map
    rad_zone = [];
    for z = 1:nzone
        nr = numel(Qnom{z});
        rad_zone = [rad_zone; z*ones(nr,1)];
    end
    Qnom_rad = zeros(numel(rad_zone),1);
    k = 0;
    for z = 1:nzone
        qz = Qnom{z};
        for i = 1:numel(qz)
            k = k + 1;
            Qnom_rad(k) = qz(i);
        end
    end
    
    P.Qnom = Qnom_rad;
    P.xmu = xmu;
    P.xsig = xsig;

    P.nx    = nx;
    P.nzone = nzone;
    P.nrad  = numel(rad_zone);

    P.vec_idx_zone = vec_idx_zone(:);
    P.rad_zone     = rad_zone;

    P.idx_x    = 1:nx;
    P.idx_Trad = nx + (1:P.nrad);
    P.idx_Tsen = nx + P.nrad + (1:nzone);
    P.idx_Tw   = nx + P.nrad + nzone + (1:P.nrad);


    % --- indices of Tsen in the FULL state vector
    idx_Tsen_global = nx + P.nrad + (1:nzone);

    % --- store physical setpoints too (recommended)
    P.Tdes_phys = Tdes(:);

    % --- zone-wise sigmas for sensor temps
    P.sig_Tsen = xsig(idx_Tsen_global);   % [nzone x 1]

    % --- normalized Tdes per radiator (matched to its zone scaling)
    Tdes_norm = zeros(P.nrad,1);
    for r = 1:P.nrad
        z = P.rad_zone(r);
        idx = idx_Tsen_global(z);
        Tdes_norm(r) = (Tdes(r) - xmu(idx)) / xsig(idx);
    end
    P.Tdes = Tdes_norm;

    P.kappa_phys = 5;
    P.sig_Tsen   = xsig(nx + P.nrad + (1:nzone));      % nzone x 1
    P.kappa_norm_zone = P.kappa_phys .* P.sig_Tsen;    % nzone x 1
    P.two_norm_zone   = 2 ./ P.sig_Tsen;               % nzone x 1   (2 K shift)

    P.e_gain  = 0.25;
    P.dT_gain = 0.10;
end

