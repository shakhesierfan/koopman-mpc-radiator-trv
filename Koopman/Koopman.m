clear; clc;

import casadi.*
load('T0.mat')
load('solar.mat')
lat = 51.652350; 
lon = 5.649608;


Ts_hrs = 0.25;

Np = 10;
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

%% Collecting Data

M = 500;
Ts = Ts_hrs*60*60;          % sampling time [s]
Tend = 2*24*3600;           % per experiment duration [s]
tgrid = (0:Ts:Tend).';
N = numel(tgrid);

X_cell  = cell(M,1);
X1_cell = cell(M,1);
U_cell  = cell(M,1);
D_cell  = cell(M,1);


x_sim = [];


for j = 1:M

    x0    = (10 + (24 - 10) * rand(nx,1))     + T_to_K;
    Trad0 = (10 + (24 - 10) * rand(n_rad, 1)) + T_to_K;
    Tw0   = (10 + (24 - 10) * rand(n_rad, 1)) + T_to_K;
    Tsen0 = (10 + (24 - 10) * rand(nzone,1))  + T_to_K;
    s0 = [x0; Trad0; Tsen0; Tw0];

    umin = 25 + T_to_K;
    umax = 95 + T_to_K;
    
    if mod(j, 3) == 0
        u = koopman_u_slowsteps_ramped(tgrid, Ts, umin, umax);
    elseif mod(j, 3) == 1
        u = generate_sin_random_peaks(tgrid, Ts, umin, umax);
    else
        u = make_piecewise_exciting_u(tgrid, Ts, umin, umax);
    end

    dk  = [T0(1:N-1).'; Int_Gain(:,1:N-1); sol(:,1:N-1)];

    sSim = simulate_full_nonlinear(tgrid, u, s0, rhs_base, T0, Int_Gain.', sol.');
        
    Xk  = sSim(:,1:end-1);
    Xk1 = sSim(:,2:end);
    Uk  = u(1:end-1).';

    X_cell{j}  = Xk;
    X1_cell{j} = Xk1;
    U_cell{j}  = Uk;
    D_cell{j}  = dk;

end

%%

rng(1);
perm = randperm(M);
Mval = round(0.2*M);
val_idx = perm(1:Mval);
tr_idx  = perm(Mval+1:end);

Xtr  = cat(2, X_cell{tr_idx});
X1tr = cat(2, X1_cell{tr_idx});
Utr  = cat(2, U_cell{tr_idx});
Dtr  = cat(2, D_cell{tr_idx});

Xva  = X_cell(val_idx);
X1va = X1_cell(val_idx);
Uva  = U_cell(val_idx);
Dva  = D_cell(val_idx);

for i = 1:Mval
    x_sim = [x_sim, [Xva{i}, X1va{i}(:, end)]];
end
%% ---------- KOOPMAN

% Normalize inputs/disturbances from training only
umu  = mean(Utr,2);  usig = std(Utr,0,2)+1e-12;
dmu  = mean(Dtr,2);  dsig = std(Dtr,0,2)+1e-12;
xmu  = mean(Xtr,2);  xsig = std(Xtr,0,2)+1e-12;

Xtrn  = (Xtr  - xmu)./xsig;
X1trn = (X1tr - xmu)./xsig;

%%

liftP = build_lift_params(nx, nzone, Qnom, vec_idx_zone, Tdes, xmu, xsig, Az, Vzone);
liftP.delta_N = delta_N;
liftP.xmu = xmu;
liftP.xsig = xsig;
liftP.umu = umu;
liftP.usig = usig;
liftP.dmu = dmu;
liftP.dsig = dsig;

Ztr  = lift_poly(Xtrn,  liftP, false);   % true = linlift
Z1tr = lift_poly(X1trn, liftP, false);


tol  = 1e-12;
zsig = std(Ztr,0,2);
zmu  = mean(Ztr,2);
badZ = zsig < tol;

zmu(1)  = 0;
zsig(1) = 1;


zsig(badZ) = 1;


Ztrs = (Ztr - zmu)./zsig;
Ztrs(badZ,:) = 0;

Z1trs = (Z1tr - zmu)./zsig;
Z1trs(badZ,:) = 0;

Ztrs(1,:)  = 1;
Z1trs(1,:) = 1;

Un = (Utr - umu)./usig;
Dn = (Dtr - dmu)./dsig;

nObs = size(Ztrs,1);

%Iu = [ 2+nx:2+nx+n_rad-1   (n_states + n_rad + 1: n_states + 2*n_rad) + 1 ];
Iu = (2:nObs); %nObs-n_rad + 1

Zu = Ztrs(Iu,:) .* Un;      % this matches rollout: z(Iu)*un


Theta_bi = [Ztrs; Un; Zu; Dn];

lambdaC = 1e-6;
Ck_bi = (Xtrn * Ztrs.') / (Ztrs*Ztrs.' + lambdaC*eye(nObs));


[xhat_best_bi, best_lam] = find_best_lam(Z1trs, ...
    Theta_bi, Iu, Ck_bi, Xva, X1va, Uva, Dva, Mval, n_states, liftP,...
    xmu, xsig, umu, usig, dmu, dsig, zmu, zsig, false, Np, badZ);


%best_lam = 3.52e-01; 4.88e-01; %1.07e-03;

Kmat_bi = fit_edmd_ridge_svd(Z1trs, Theta_bi, best_lam);

idx = 0;
nObs = size(Ztrs,1);
A_bi = real(Kmat_bi(:, idx+1:idx+nObs)); idx = idx+nObs;
B_bi = real(Kmat_bi(:, idx+1)); idx = idx+1;
nZu  = numel(Iu);
N_bi = real(Kmat_bi(:, idx+1:idx+nZu)); idx = idx+nZu;

BT0_bi  = real(Kmat_bi(:, idx+1)); idx = idx+1;
nInt    = liftP.nzone;
Bint_bi = real(Kmat_bi(:, idx+1:idx+nInt)); idx = idx+nInt;
Bsol_bi = real(Kmat_bi(:, idx+1:end));

% enforce constant observable
A_bi(1,:)=0; A_bi(1,1)=1;
B_bi(1)=0; N_bi(1,:)=0;
BT0_bi(1)=0; Bint_bi(1,:)=0; Bsol_bi(1,:)=0;






%% Koopman (Linear)


% ===== Linear Koopman (nonlinear lifting) + ridge scan, using rollout_rmse_bilinear =====

% Lift (nonlinear) on TRAINING data
Ztr_lin  = lift_poly(Xtrn,  liftP, false);   % false = nonlinear lifting
Z1tr_lin = lift_poly(X1trn, liftP, false);

% Scale lifted observables (TRAINING)
tol = 1e-8;

zsig_lin = std(Ztr_lin,0,2);
zmu_lin = mean(Ztr_lin,2);
badZ_lin = zsig_lin < tol;

zmu_lin(1)  = 0;
zsig_lin(1) = 1;


zsig_lin(badZ_lin) = 1;          % avoid blow-up
Ztrs_lin = (Ztr_lin - zmu_lin)./zsig_lin;
Ztrs_lin(badZ_lin,:) = 0;        % keep dimension, kill feature

Ztrs_lin(1,:)  = 1;

Z1trs_lin = (Z1tr_lin - zmu_lin)./zsig_lin;
Z1trs_lin(badZ_lin,:) = 0;
Z1trs_lin(1,:) = 1;


Un_tr = (Utr - umu)./usig;
Dn_tr = (Dtr - dmu)./dsig;

Theta_lin_tr = [Ztrs_lin; Un_tr; Dn_tr];

% sizes
nObs_lin = size(Ztrs_lin,1);

lambdaC = 1e-6;
C_lin_K = (Xtrn * Ztrs_lin.') / (Ztrs_lin*Ztrs_lin.' + lambdaC*eye(nObs_lin));

[xhat_best_lin, best_lam_lin] = find_best_lam(Z1trs_lin, ...
    Theta_lin_tr, [], C_lin_K, Xva, X1va, Uva, Dva, Mval, n_states, liftP,...
    xmu, xsig, umu, usig, dmu, dsig, zmu_lin, zsig_lin, false, Np, badZ_lin);


Kmat_best_lin = fit_edmd_ridge_svd(Z1trs_lin, Theta_lin_tr, best_lam_lin);

idx = 0;
nInt = liftP.nzone;

A_lin_K  = real(Kmat_best_lin(:, idx+1:idx+nObs_lin)); idx = idx+nObs_lin;
B_lin_K  = real(Kmat_best_lin(:, idx+1));             idx = idx+1;

BT0_lin_K  = real(Kmat_best_lin(:, idx+1)); idx = idx+1;
Bint_lin_K = real(Kmat_best_lin(:, idx+1:idx+nInt)); idx = idx+nInt;
Bsol_lin_K = real(Kmat_best_lin(:, idx+1:end));

A_lin_K(1,:)=0;  A_lin_K(1,1)=1;
B_lin_K(1)=0;
BT0_lin_K(1)=0;
Bint_lin_K(1,:)=0;
Bsol_lin_K(1,:)=0;



%% Koopman (Linear lifting)


% ----- linlift training lifted data -----
Ztr_linlift  = lift_poly(Xtrn,  liftP, true);
Z1tr_linlift = lift_poly(X1trn, liftP, true);

tol = 1e-8;

zsig_linlift = std(Ztr_linlift,0,2);
badZ_linlift = zsig_linlift < tol;

zsig_linlift(badZ_linlift) = 1;
zmu_linlift = mean(Ztr_linlift,2);

zmu_linlift(1)  = 0;
zsig_linlift(1) = 1;


Ztrs_linlift = (Ztr_linlift - zmu_linlift)./zsig_linlift;
Ztrs_linlift(badZ_linlift,:) = 0;
Ztrs_linlift(1,:) = 1;

Z1trs_linlift = (Z1tr_linlift - zmu_linlift)./zsig_linlift;
Z1trs_linlift(badZ_linlift,:) = 0;
Z1trs_linlift(1,:) = 1;


Un_tr = (Utr - umu)./usig;
Dn_tr = (Dtr - dmu)./dsig;

Theta_linlift_tr = [Ztrs_linlift; Un_tr; Dn_tr];

nObs_linlift = size(Ztrs_linlift,1);
lambdaC = 1e-6;
C_linlift_K = (Xtrn * Ztrs_linlift.') / (Ztrs_linlift*Ztrs_linlift.' + lambdaC*eye(nObs_linlift));


[xhat_best_linlift, best_lam_linlift] = find_best_lam(Z1trs_linlift, ...
    Theta_linlift_tr, [], C_linlift_K, Xva, X1va, Uva, Dva, Mval, n_states, liftP,...
    xmu, xsig, umu, usig, dmu, dsig, zmu_linlift, zsig_linlift, true, Np, badZ_linlift);


Kmat_best_linlift = fit_edmd_ridge_svd(Z1trs_linlift, Theta_linlift_tr, best_lam_linlift);

idx = 0;
nInt = liftP.nzone;

A_linlift_K  = real(Kmat_best_linlift(:, idx+1:idx+nObs_linlift)); idx = idx+nObs_linlift;
B_linlift_K  = real(Kmat_best_linlift(:, idx+1));                  idx = idx+1;

BT0_linlift_K  = real(Kmat_best_linlift(:, idx+1)); idx = idx+1;
Bint_linlift_K = real(Kmat_best_linlift(:, idx+1:idx+nInt)); idx = idx+nInt;
Bsol_linlift_K = real(Kmat_best_linlift(:, idx+1:end));

% constant observable
A_linlift_K(1,:)=0;  A_linlift_K(1,1)=1;
B_linlift_K(1)=0;
BT0_linlift_K(1)=0;
Bint_linlift_K(1,:)=0;
Bsol_linlift_K(1,:)=0;






%%  COMPARISON



Lin = linearize_ct_analytical(A, Brad, BT0, Bint, Bsol, vec_idx_zone, Az, Qnom, delta_N, Tdes, Tdel_sec);

%{
t_T0  = (0:size(T0, 1)-1) * (Ts_hrs*3600); 
t_Int = (0:size(Int_Gain, 2)-1) * (Ts_hrs*3600); 
t_sol = (0:size(sol, 2)-1) * (Ts_hrs*3600); 


T0_fun  = @(t) interp1(t_T0,  T0,       t, 'previous', 'extrap');
Int_fun = @(t) interp1(t_Int, Int_Gain.',t, 'previous', 'extrap').';
Sol_fun = @(t) interp1(t_sol, sol.',     t, 'previous', 'extrap').';

%}


%s_op_k  = (25 + T_to_K)*ones(n_states,1);
%u_op_k = 50 + T_to_K;
%T0_op_k = 10 + T_to_K;
%Int_op_k = 5*ones(nzone, 1);
%Sol_op_k = 20*ones(size(Bsol,2), 1);


A_ct    = Lin.A(xmu, umu, dmu(1), dmu(2:nzone+1), dmu(nzone+2:end));
B_ct    = Lin.B(xmu, umu, dmu(1), dmu(2:nzone+1), dmu(nzone+2:end));
ET0_ct  = Lin.ET0(xmu, umu, dmu(1), dmu(2:nzone+1), dmu(nzone+2:end));
EInt_ct = Lin.EInt(xmu, umu, dmu(1), dmu(2:nzone+1), dmu(nzone+2:end));
ESol_ct = Lin.ESol(xmu, umu, dmu(1), dmu(2:nzone+1), dmu(nzone+2:end));


eps_stab = 1e-3;
alpha = max(real(eig(A_ct)));

if alpha > -eps_stab
    A_ct = A_ct - (alpha + eps_stab)*eye(size(A_ct));
end

nu   = size(B_ct,2);
iT0  = size(ET0_ct,2);
nInt = size(EInt_ct,2);
nSol = size(ESol_ct,2);
nE   = iT0 + nInt + nSol;




fk = Lin.f(xmu, umu, dmu(1), dmu(2:nzone+1), dmu(nzone+2:end));
ck = fk - A_ct*xmu - B_ct*umu - ET0_ct*dmu(1) - EInt_ct*dmu(2:nzone+1) - ESol_ct*dmu(nzone+2:end);

B_st_ct = [B_ct ET0_ct EInt_ct ESol_ct ck];


M_exp = expm([A_ct B_st_ct;
          zeros(nu+nE+1, n_states+nu+nE+1)] * Ts);

Ad_init = M_exp(1:n_states, 1:n_states);
idx = n_states + 1;
Bd_init = M_exp(1:n_states, idx:idx+nu-1); idx = idx + nu;
ET0d_init = M_exp(1:n_states, idx:idx+iT0-1); idx = idx + iT0;
EIntd_init = M_exp(1:n_states, idx:idx+nInt-1); idx = idx + nInt;
ESold_init= M_exp(1:n_states, idx:idx+nSol-1); idx = idx + nSol;
cd_init   = M_exp(1:n_states, idx);


%%
%X = zeros(n_states, N*n_comp);
%X_init = zeros(n_states, N*n_comp);


nExp = numel(Xva);
se_LTV = 0; 
se_LTI = 0;
cnt = 0;

xhat_sim_LTI = zeros(n_states, Mval*(size(Xva{1},2) + 1));
xhat_sim_LTV = zeros(n_states, Mval*(size(Xva{1},2) + 1));

idx_sim = 0;
for e = 1:nExp
    Xk  = Xva{e};
    Xk1 = X1va{e};
    U = Uva{e};
    D = Dva{e};

    %x0 = Xk(:,1);

    Nsteps = min([size(U,2), size(D,2), size(Xk,2), size(Xk1,2)]);

    
    idx_sim = idx_sim + 1;
    %xhat_sim_LTI(:, idx_sim) = Xk(:,1);
    %xhat_sim_LTV(:, idx_sim) = Xk(:,1);
    
    for k = 1:Nsteps

        maxH = min(Np, size(Xk,2) - k);

        if maxH <= 0
            %max = 1;
        %elseif maxH < 0
            continue;
        end
        %    continue;
        %end

        xhat_LTV = Xk(:,k);
        xhat_LTI = Xk(:,k);
        for h = 1:maxH
            kk = k + h - 1;

    
            T0k  = D(1 ,kk);
            Intk = D(2:1+nzone ,kk);
            Solk = D(2+nzone:end ,kk);
        
            s_op_k = xhat_LTV;
            u_op_k = U(kk);
    
            T0_op_k  = T0k;
            Int_op_k = Intk;
            Sol_op_k = Solk;
        
            A_ct_LTV    = Lin.A(s_op_k, u_op_k, T0_op_k, Int_op_k, Sol_op_k);
            B_ct_LTV    = Lin.B(s_op_k, u_op_k, T0_op_k, Int_op_k, Sol_op_k);
            ET0_ct_LTV  = Lin.ET0(s_op_k, u_op_k, T0_op_k, Int_op_k, Sol_op_k);
            EInt_ct_LTV = Lin.EInt(s_op_k, u_op_k, T0_op_k, Int_op_k, Sol_op_k);
            ESol_ct_LTV = Lin.ESol(s_op_k, u_op_k, T0_op_k, Int_op_k, Sol_op_k);
    
    
            eps_stab = 1e-3;
            alpha = max(real(eig(A_ct_LTV)));
            
            if alpha > -eps_stab
                A_ct_LTV = A_ct_LTV - (alpha + eps_stab)*eye(size(A_ct_LTV));
            end
    
            
    
            fk_LTV = Lin.f(s_op_k, u_op_k, T0_op_k, Int_op_k, Sol_op_k);
            ck_LTV = fk_LTV - A_ct_LTV*s_op_k - B_ct_LTV*u_op_k - ET0_ct_LTV*T0_op_k ...
                - EInt_ct_LTV*Int_op_k - ESol_ct_LTV*Sol_op_k;
    
    
            B_st = [B_ct_LTV ET0_ct_LTV EInt_ct_LTV ESol_ct_LTV ck_LTV];
            
            M_exp = expm([A_ct_LTV B_st;
              zeros(nu+nE+1, n_states+nu+nE+1)] * Ts);
                
            Ad = M_exp(1:n_states, 1:n_states);
            idx = n_states + 1;
            Bd = M_exp(1:n_states, idx:idx+nu-1); idx = idx + nu;
            ET0d = M_exp(1:n_states, idx:idx+iT0-1); idx = idx + iT0;
            EIntd= M_exp(1:n_states, idx:idx+nInt-1); idx = idx + nInt;
            ESold= M_exp(1:n_states, idx:idx+nSol-1); idx = idx + nSol;
            cd   = M_exp(1:n_states, idx);
        
            xhat_LTV = Ad*xhat_LTV+ Bd*u_op_k + ET0d*T0k + EIntd*Intk + ESold*Solk + cd;
            xhat_LTI = Ad_init*xhat_LTI + Bd_init*u_op_k + ET0d_init*T0k + EIntd_init*Intk + ESold_init*Solk + cd_init;
        
    
    
    
            xtrue = Xk1(:, kk);

            err_LTV = xhat_LTV - xtrue;
            err_LTI = xhat_LTI - xtrue;
            se_LTV  = se_LTV + sum(err_LTV(:).^2);
            se_LTI  = se_LTI + sum(err_LTI(:).^2);

            cnt = cnt + numel(err_LTV);
        end



    end

end
rmse_LTV = sqrt(se_LTV / max(cnt,1));
rmse_LTI = sqrt(se_LTI / max(cnt,1));


fprintf('\n[LTV Linear]\n');
fprintf('  RMSE global        : %.6f\n', rmse_LTV);

fprintf('\n[LTI Linear]\n');
fprintf('  RMSE global        : %.6f\n', rmse_LTI);


%% ----- Plotting ---------

N = 1000;
figure;
plot(x_sim(end,N:2*N), 'LineWidth',1.5); hold on;
plot(xhat_best_bi(end,N:2*N), 'LineWidth',1.5); hold on;
plot(xhat_best_lin(end,N:2*N), 'LineWidth',1.5); hold on;
plot(xhat_best_linlift(end,N:2*N), '--', 'LineWidth',1.5); hold on;
%plot(xhat_sim_LTI(1,1:N), '--', 'LineWidth',1.5, 'Color','g'); hold on;
plot(xhat_sim_LTV(end,N:2*N), '*', 'LineWidth',1.5); hold on;
legend('Nonlinear','Koopman', 'Koopman (linear)', 'Koopman (linear lifting)', 'LTV Lin');




%{
figure;
plot(tgrid/3600, sNLs(1,1:N), 'LineWidth',1.5); hold on;
plot(tgrid/3600, sHat(1,1:N), '--', 'LineWidth',1.5); hold on;
plot(tgrid/3600, sHat_lin(1,1:N), '--', 'LineWidth',1.5, 'Color','g'); hold on;
plot(tgrid/3600, X(1,1:N), '*', 'LineWidth',1.5); hold on;
plot(tgrid/3600, X_init(1,1:N), '*', 'LineWidth',1.5);

grid on; xlabel('Time [h]'); ylabel('State component 1');
legend('Nonlinear','Koopman', 'Koopman (linear)', 'LTV Lin', 'LTI Lin');

figure;
plot(tgrid/3600, sNLs(42,1:N), 'LineWidth',1.5); hold on;
plot(tgrid/3600, sHat(42,1:N), 'o', 'LineWidth',1.5); hold on;
plot(tgrid/3600, sHat_lin(42,1:N), 'o', 'LineWidth',1.5, 'Color','g'); hold on;
plot(tgrid/3600, X(42,1:N), '*', 'LineWidth',1.5); hold on;
plot(tgrid/3600, X_init(42,1:N), '*', 'LineWidth',1.5);


grid on; xlabel('Time [h]'); ylabel('State component 1');
legend('Nonlinear','Koopman', 'Koopman (linear)', 'LTV Lin', 'LTI Lin');

%}



%% == FUNCTIONS ==


function sSim = simulate_full_nonlinear(tgrid, u, s0, rhs_full, T0, Int, Sol)
% Fully ZOH simulation:
%   u, T0, IntG, Sol are held constant on each [kTs,(k+1)Ts] interval.

    opts = odeset('RelTol',1e-7,'AbsTol',1e-12);

    sSim = s0;
    Ts = tgrid(2) - tgrid(1);

    for k = 1:(numel(tgrid)-1)

        % --- Freeze inputs at the LEFT endpoint (ZOH) ---
        %tk   = tgrid(k);
        uk   = u(k);

        T0k   = T0(k);
        Intk  = Int(k,:).';
        Solk  = Sol(k,:).';

        % Constant functions over this interval
        Tsup_fun_k = @(t) uk;
        T0_fun_k  = @(t) T0k;
        Int_fun_k = @(t) Intk;
        Sol_fun_k = @(t) Solk;

        % Build an RHS that uses constants (no time-varying disturbance calls)
        rhs_k = @(t,s) rhs_full(0 , s, Tsup_fun_k, T0_fun_k, Int_fun_k, Sol_fun_k);

        % Integrate one step
        [~, sTemp] = ode45(rhs_k, [0 Ts], s0, opts);  % local time 0..Ts
        s0 = sTemp(end,:).';
        sSim = [sSim s0];
    end
    
    if ~isreal(sSim)
        error('simulate_full_nonlinear: complex states detected');
    end
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
        qz = Qnom{z};          % e.g. [1700 1700 ...] or [] 
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

    % --- store physical setpoints
    P.Tdes_phys = Tdes(:);

    % --- zone-wise sigmas for sensor temps
    P.sig_Tsen = xsig(idx_Tsen_global);   % [nzone x 1]

    % --- normalized Tdes per radiator
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













%% Linearized 

function Lin = linearize_ct_analytical(A, Brad, BT0, Bint, Bsol, vec_idx_zone, Az, Qnom, delta_N, Tdes, Tdel_sec)
    nx    = size(A, 1);
    nzone = numel(vec_idx_zone);
    n_rad = numel([Qnom{:}]);
    nInt = nzone;
    nSol = size(Bsol,2);
    
    
    %t  = sym('t','real');
    s  = sym('s',[nx+2*n_rad+nzone, 1],'real');
    u  = sym('u','real');
    T0 = sym('T0','real');
    IntG = sym('IntG', [nInt,1],'real');
    Sol  = sym('Sol', [nSol,1],'real');
    
    x    = s(1:nx);
    Trad = s(nx+1:nx+n_rad);
    Tsen = s(nx+n_rad+1:nx+n_rad+nzone);
    Tw   = s(nx+n_rad+nzone+1:end);
    
    
    q_rad = compute_q_rad_ct(Qnom, Trad, x, Az, delta_N, vec_idx_zone, u);
    
    
    [dTrad, dTsen, dTw] = compute_T_rad_ct(x, vec_idx_zone, Tdes, u, Qnom, Trad, Tsen, Tdel_sec, delta_N, Tw);
    
    
    dx = A*x + BT0*T0 + Bint*IntG + Bsol*Sol + Brad*q_rad;
    ds = [dx; dTrad; dTsen; dTw];
    
    A_jac    = jacobian(ds, s);
    B_jac    = jacobian(ds, u);
    ET0_jac  = jacobian(ds, T0);
    EInt_jac = jacobian(ds, IntG);
    ESol_jac = jacobian(ds, Sol);
    
    vars = {s,u,T0,IntG,Sol};
    
    Lin.f    = matlabFunction(ds,      'Vars', vars);
    Lin.A    = matlabFunction(A_jac,   'Vars', vars);
    Lin.B    = matlabFunction(B_jac,   'Vars', vars);
    Lin.ET0  = matlabFunction(ET0_jac, 'Vars', vars);
    Lin.EInt = matlabFunction(EInt_jac,'Vars', vars);
    Lin.ESol = matlabFunction(ESol_jac,'Vars', vars);

end

%{
function delta = lmtd_stable(d1, d2)
    epsd  = sym(1e-6);
    epsln = sym(1e-12);

    d1p = sqrt(d1.^2 + epsd^2);
    d2p = sqrt(d2.^2 + epsd^2);

    r = d1p./ d2p;
    y = r - 1;
    den = log(1 + y);

    phi = y ./ (den + epsln);

    delta = d2p .* phi;
end
%}


function u = koopman_u_slowsteps_zoh(tgrid, Ts, umin, umax)

    tgrid = tgrid(:);
    t0 = tgrid(1);
    tf = tgrid(end);

    N = floor((tf - t0)/Ts) + 1;     % number of Ts samples
    n = (0:N-1)';

    M = max(6, min(24, floor(N/20)));   % adaptive default; change if you know better

    K = ceil(N / M);

    levels = 2*rand(K,1) - 1;

    mid  = (umin + umax)/2;
    half = (umax - umin)/2;
    levels = mid + half*levels;

    uS = repelem(levels, M);
    uS = uS(1:N);

    idx = floor((tgrid - t0)/Ts) + 1;
    idx = min(max(idx,1),N);
    u = uS(idx);

    u = reshape(u, size(tgrid));
end


function [u, k] = generate_sin_random_peaks(tgrid, Ts, umin, umax, varargin)

    kmin = 1;
    kmax = 5;
    seed = [];

    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'mincycles'
                kmin = varargin{i+1};
            case 'maxcycles'
                kmax = varargin{i+1};
            case 'seed'
                seed = varargin{i+1};
            otherwise
                error('Unknown option: %s', varargin{i});
        end
    end

    if ~isempty(seed)
        rng(seed);
    end

    k = randi([kmin, kmax]);

    umid = (umax + umin)/2;
    A    = (umax - umin)/2;


    T_end = tgrid(end) - tgrid(1);
    omega = 2*pi*k / T_end;

    t0 = tgrid(1);
    phi = -pi/2 - omega*t0;

    u = umid + A*sin(omega*tgrid + phi);

    u = min(max(u, umin), umax);
end

function u = koopman_u_slowsteps_ramped(tgrid, Ts, umin, umax, M, ramp_frac)

    if nargin < 5 || isempty(M)
        M = 24;
    end
    if nargin < 6 || isempty(ramp_frac)
        ramp_frac = 0.2;
    end

    tgrid = tgrid(:);
    t0 = tgrid(1);
    tf = tgrid(end);

    N = floor((tf - t0)/Ts) + 1;
    K = ceil(N / M);

    % random plateau levels in [umin, umax]
    mid  = (umin + umax)/2;
    half = (umax - umin)/2;
    levels = mid + half*(2*rand(K,1) - 1);

    start_frac = 0.05;
    levels(1) = umin + start_frac*(umax - umin);

    ramp_len = max(1, min(M-1, round(ramp_frac*M)));

    uS = zeros(N,1);
    k0 = 1;

    for k = 1:K
        a = levels(k);
        if k < K
            b = levels(k+1);
        else
            b = levels(k);
        end

        k1 = min(N, k0 + M - 1);
        L  = k1 - k0 + 1;

        r = min(ramp_len, max(0, L-1));
        hold_len = L - r;

        uS(k0 : k0+hold_len-1) = a;

        if r > 0
            s = (1:r)'/r;
            uS(k0+hold_len : k1) = (1 - s).*a + s.*b;
        end

        k0 = k1 + 1;
        if k0 > N, break; end
    end

    idx = floor((tgrid - t0)/Ts) + 1;
    idx = min(max(idx,1),N);
    u = uS(idx);
    u = reshape(u, size(tgrid));
end

function u = make_piecewise_exciting_u(tgrid, Ts, umin, umax)

    if umax <= umin
        error('umax must be > umin.');
    end
    tgrid = tgrid(:);
    N = numel(tgrid);
    if N == 0
        u = zeros(0,1);
        return;
    end

    kidx = round(tgrid./Ts);           % sample index for each tgrid entry
    kmin = kidx(1);
    kmax = kidx(end);
    K = kmax - kmin + 1;               % number of sample instants covered

    u_k = zeros(K,1);

    slow_min_h = 2.0;
    slow_max_h = 6.0;

    fast_min_h = 0.25;
    fast_max_h = 0.50;

    fast_amp = 0.08*(umax-umin);
    du_max   = inf; 
    % ------------------------------------------------------

    slow_min = max(1, round(slow_min_h*3600/Ts));
    slow_max = max(slow_min, round(slow_max_h*3600/Ts));
    fast_min = max(1, round(fast_min_h*3600/Ts));
    fast_max = max(fast_min, round(fast_max_h*3600/Ts));

    u_k(1) = umin + (umax-umin)*rand;

    k = 1;
    while k <= K
        slow_dwell = randi([slow_min, slow_max]);
        u_slow = umin + (umax-umin)*rand;

        kend = min(k + slow_dwell - 1, K);

        kk = k;
        while kk <= kend
            fast_dwell = randi([fast_min, fast_max]);
            fast_end = min(kk + fast_dwell - 1, kend);

            dither = fast_amp*(2*rand - 1);      % uniform in [-amp, +amp]
            level  = u_slow + dither;

            u_k(kk:fast_end) = level;            % constant over this block
            kk = fast_end + 1;
        end

        k = kend + 1;
    end

    u_k = min(max(u_k, umin), umax);

    if isfinite(du_max)
        for k = 2:K
            u_k(k) = min(max(u_k(k), u_k(k-1) - du_max), u_k(k-1) + du_max);
        end
    end

    u = zeros(N,1);
    for i = 1:N
        kk = kidx(i) - kmin + 1;   % map to 1..K
        kk = min(max(kk, 1), K);
        u(i) = u_k(kk);
    end
end

function Kmat = fit_edmd_ridge_svd(Z1, Theta, lam)

    [U,S,V] = svd(Theta,'econ');
    s = diag(S);

    % Truncate tiny singular values (relative)
    tol = 1e-12 * max(s);
    keep = s > tol;

    U = U(:,keep); V = V(:,keep); s = s(keep);

    Theta_pinv = V * diag(s ./ (s.^2 + lam)) * U';
    Kmat = Z1 * Theta_pinv;

end


function [rmse, xhat_sim] = rollout_rmse_bilinear_Np( ...
    Ak,Bk,Nk,Iu,ET0k,EIntk,ESolk,Ck, ...
    Xva,X1va,Uva,Dva, xmu,xsig, umu,usig, dmu,dsig, ...
    liftP, zmu,zsig, useLinearLift, Np, badZ)

    if nargin < 22 || isempty(Np)
        Np = 1;
    end

    nExp = numel(Xva);
    se = 0; cnt = 0;

    xhat_sim = [];
    idx_sim = 0;

    for e = 1:nExp
        Xk  = Xva{e};
        Xk1 = X1va{e};
        U   = Uva{e};
        D   = Dva{e};

        idx_sim = idx_sim + 1;
        xhat_sim(:, idx_sim) = Xk(:,1);

        Nsteps = min([size(U,2), size(D,2), size(Xk,2), size(Xk1,2)]);
        if Nsteps <= 0, continue; end


        for k = 1:Nsteps
            xk0n = (Xk(:,k) - xmu)./xsig;
            z0   = lift_poly(xk0n, liftP, useLinearLift);
            zk   = (z0 - zmu)./zsig;
            zk(badZ) = 0;
            zk(1) = 1;

            un = (U(k) - umu)./usig;
            dn = (D(:,k) - dmu)./dsig;

            T0n  = dn(1);
            Intn = dn(2:1+liftP.nzone);
            Soln = dn(2+liftP.nzone:end);

            if ~isempty(Iu)
                zu = zk(Iu) * un;
                zk_seq = Ak*zk + Bk*un + Nk*zu + ET0k*T0n + EIntk*Intn + ESolk*Soln;
            else
                zk_seq = Ak*zk + Bk*un + ET0k*T0n + EIntk*Intn + ESolk*Soln;
            end
            zk_seq(badZ) = 0;
            zk_seq(1) = 1;

            xhatn = Ck*zk_seq;
            xhat  = xhatn.*xsig + xmu;

            idx_sim = idx_sim + 1;

            xhat_sim(:, idx_sim) = xhat;
            

        end

        Kmax = Nsteps;  % last "current" index
        for k0 = 1:Kmax
            maxH = min(Np, size(Xk,2) - k0);

            if maxH <= 0
                %max = 1;
            %elseif maxH < 0
                continue;
            end
            %    continue;
            %end

            xk0n = (Xk(:,k0) - xmu)./xsig;
            z0   = lift_poly(xk0n, liftP, useLinearLift);
            zk   = (z0 - zmu)./zsig;
            zk(badZ) = 0;
            zk(1) = 1;
            for h = 1:maxH
                kk = k0 + h - 1;

                un = (U(kk) - umu)./usig;
                dn = (D(:,kk) - dmu)./dsig;

                T0n  = dn(1);
                Intn = dn(2:1+liftP.nzone);
                Soln = dn(2+liftP.nzone:end);

                if ~isempty(Iu)
                    zu = zk(Iu) * un;
                    zk = Ak*zk + Bk*un + Nk*zu + ET0k*T0n + EIntk*Intn + ESolk*Soln;
                else
                    zk = Ak*zk + Bk*un + ET0k*T0n + EIntk*Intn + ESolk*Soln;
                end
                zk(badZ) = 0;

                zk(1) = 1;

                xhatn = Ck*zk;
                xhat  = xhatn.*xsig + xmu;

                % true future state at time k0+h
                xtrue = Xk1(:, kk); %~~~~~~~~~~~~~~~~~~~~~

                err = xhat - xtrue;
                se  = se + sum(err(:).^2);
                cnt = cnt + numel(err);
            end
        end

    end

    rmse = sqrt(se / max(cnt,1));
    xhat_sim = xhat_sim(:);
end



function [xhat_best_bi, best_lam] = find_best_lam(Z1trs, Theta, Iu, Ck, Xva, X1va,Uva,Dva, Mval, n_states, liftP,...
    xmu, xsig, umu, usig, dmu, dsig, zmu, zsig, useLinearLift, Np, badZ)

    smax = svds(Theta,1);
    lams = 2.56e+00; %logspace(-12, -1, 40) * (smax^2);
    %lams = 3.52e-01;
    rmse_val = zeros(size(lams));
    
    %xhat_sim = zeros(Mval*(size(Xva{1},2) + 1)*n_states, 1);
    best_rmse_temp = 100;
    
    for i=1:numel(lams)
        lam = lams(i);
    
        % fit bilinear EDMD
        Kmat = fit_edmd_ridge_svd(Z1trs, Theta, lam);
    
        idx = 0;
    
        nObs = size(Z1trs,1);
        A = real(Kmat(:, idx+1:idx+nObs)); idx = idx+nObs;
    
        B = real(Kmat(:, idx+1)); idx = idx+1;
    
        nZu = numel(Iu);
        if nZu~= 0
            N = real(Kmat(:, idx+1:idx+nZu)); idx = idx+nZu;
        else
            N = 0;
        end
    
        % disturbance blocks
        BT0  = real(Kmat(:, idx+1)); idx = idx+1;
        nInt    = liftP.nzone;
        Bint = real(Kmat(:, idx+1:idx+nInt)); idx = idx+nInt;
        Bsol = real(Kmat(:, idx+1:end));
    
        % enforce constant observable
        A(1,:)=0; A(1,1)=1;
        B(1)=0; N(1,:)=0;
        BT0(1)=0; Bint(1,:)=0; Bsol(1,:)=0;
        
        % rollout bilinear RMSE on validation
        [rmse_val(i), xhat_sim_temp] = rollout_rmse_bilinear_Np( ...
            A,B,N,Iu,BT0,Bint,Bsol,Ck, ...
            Xva,X1va,Uva,Dva, xmu,xsig, umu,usig, dmu,dsig, ...
            liftP, zmu,zsig, useLinearLift, Np, badZ);
        if rmse_val(i) <= best_rmse_temp
            best_rmse_temp = rmse_val(i);
            xhat_sim = xhat_sim_temp;
        end

        fprintf("lam=%9.2e  val_RMSE=%g\n", lam, rmse_val(i));
    end
    
    [best_rmse, best_i] = min(rmse_val);
    best_lam = lams(best_i);
    xhat_best_bi = reshape( ...
        xhat_sim, ...
        n_states, ...
        Mval * (size(Xva{1}, 2) + 1) ...
    );
    
    fprintf("BEST lam (nonlinear) = %.3e, val_RMSE=%g\n", best_lam, best_rmse);

end
