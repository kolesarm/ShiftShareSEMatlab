%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replication of ADH (2013)
% Input: data_input_empADH, endog_var_czADH
% Adao, Kolesar, Morales - 08/06/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Data Input
load data_input_empADH;
load endog_var_cz_all;

%% Preliminaries
% Numerical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 2; %n of period
I = 722*2; % n of regions
S = 396*2; % n of sectors
alpha = .05; %significance level of hypothesis test

%Data: controls and sample selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight = timepwt48;     %observation weights

controls = [ones(size(t2)), t2, l_shind_manuf_cbp, reg_encen, reg_escen, reg_midatl, reg_mount, reg_pacif, reg_satl, reg_wncen, reg_wscen, l_sh_popedu_c, l_sh_popfborn, l_sh_empl_f, l_sh_routine33, l_task_outsource];
[obs K] = size(controls);
K = K+1;

%Sector employment shares
share_emp_ind1 = reshape(share_emp_ind1980, S/T , I/T)';
share_emp_ind2 = reshape(share_emp_ind1990, S/T , I/T)';
sec_vec = reshape(sic87dd, S/T , I/T);

share_emp_ind = [share_emp_ind1, zeros(size(share_emp_ind2)); zeros(size(share_emp_ind1)), share_emp_ind2];
sec_vec = [sec_vec(:,1); sec_vec(:,1)];
sec_vec3d = floor(sec_vec/10);


%Adjust share matrix: Drop linearlly colinear columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ln = share_emp_ind;
B = ln'*ln;
tol = 1e-10;
if ~nnz(B) %X has no non-zeros and hence no independent columns
    shock_ind=[];
    return
end
[Q, R, E] = qr(B,0);
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = R(1);
end
%Rank estimation
r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
shock_ind=sort(E(1:r));
clearvars Q R E B diagr

share_emp_ind = ln(:, shock_ind);
sec_vec = sec_vec(shock_ind');
sec_vec3d = sec_vec3d(shock_ind');

%% Estimation

% First-Stage
[ hat_beta(1,1), SE(1,1), pvalue(1,1), CIl(1,1), CIu(1,1), CIt(1,1) ] = reg_ss( d_tradeusch_pw, d_tradeotch_pw_lag, controls, share_emp_ind, weight, sec_vec3d, alpha, 1, [] );
[ hat_beta(2,1), SE(2,1), pvalue(2,1), CIl(2,1), CIu(2,1), CIt(2,1) ] = reg_ss( d_tradeusch_pw, d_tradeotch_pw_lag, controls, share_emp_ind, weight, sec_vec3d, alpha, 0, [] );

%Reduced-Form
[ hat_beta(1,2), SE(1,2), pvalue(1,2), CIl(1,2), CIu(1,2), CIt(1,2) ] = reg_ss( d_sh_empl_mfg, d_tradeotch_pw_lag, controls, share_emp_ind, weight, sec_vec3d, alpha, 1, [] );
[ hat_beta(2,2), SE(2,2), pvalue(2,2), CIl(2,2), CIu(2,2), CIt(2,2) ] = reg_ss( d_sh_empl_mfg, d_tradeotch_pw_lag, controls, share_emp_ind, weight, sec_vec3d, alpha, 0, [] );

%2SLS
[ hat_beta(1,3), SE(1,3), pvalue(1,3), CIl(1,3), CIu(1,3), CIt(1,3) ] = ivreg_ss( d_sh_empl_mfg, d_tradeusch_pw, d_tradeotch_pw_lag, controls, share_emp_ind, weight, sec_vec3d, alpha, 1, [] );
[ hat_beta(2,3), SE(2,3), pvalue(2,3), CIl(2,3), CIu(2,3), CIt(2,3) ] = ivreg_ss( d_sh_empl_mfg, d_tradeusch_pw, d_tradeotch_pw_lag, controls, share_emp_ind, weight, sec_vec3d, alpha, 0, [] );
