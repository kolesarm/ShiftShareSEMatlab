function [hat_beta, SE, pvalue, CIl, CIu, CItype] = ivreg_ss(Yn, Xn, Zn, controls, ln, weight, sec_cluster_vec, alpha, AKMtype, beta0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AKM Inference - 2SLS
%Adao, Kolesar, Morales - 07/23/2018
%Updated on 08/01/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Description of arguments
%Yn: dependent variable
%Xn: endogenous regressor
%Zn: shift-share IV
%controls: Control matrix -- vector of ones if empty
%ln: matrix of shares used in shift-share regressor
%weight: observation weights
%sec_cluster_vec: vector of clusters -- no clustering if empty
%alpha: significance level for confidence interval
%AKMtype: 1 for AKM and 0 for AKM0

%Description of Output
% hat_beta: estimated coefficient on endogenous regressor
% SE: length of CI
% pvalue: p-value of H0: beta = beta0
% CIl: lower bound of CI
% CIu: upper bound of CI
% CI type:
%   0 - AKM,
%   1 - standard AKM0
%   2 - nonstandard AKM0 = [-Inf,CIl]U[CIu,Inf]
%   3 - nonstandard AKM0 = [-Inf,Inf]

%% Preliminearies
global S K obs

critical = norminv(1 - alpha/2,0,1);

if isempty(beta0) == 1
    beta0 = 0;
end

if isempty(controls) == 1
    controls = ones(size(Yn));
end
[obs K] = size(controls);
K = K+1;

% Adjust share matrix: Drop linearlly colinear columns
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

ln_check = ln(:, shock_ind);
colinear = (size(ln_check) == size(ln));
if sum(colinear) < 2
    disp('Share matrix has colinear columns')
    hat_beta =0; SE =0; pvalue=0; CIl=0; CIu=0; CItype=0;
    return
end

%Define variables for estimation
[obs S] = size(ln);
Mn = [controls, Zn].*( repmat(sqrt(weight),1,K) );   %Matrix of IVs
Gn = [controls, Xn].*( repmat(sqrt(weight),1,K) );   %Matrix of regressors
ln = ln.*( repmat(sqrt(weight),1,S) );               %Matrix of Shares
tildeYn = Yn.*(sqrt(weight));                        %Dependent Variable


%% OLS estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1 = Mn'*Gn;
P2 = (Mn'*Mn)^(-1);
P3 = Mn'*tildeYn;
hat_theta = ((P1'*P2*P1)^(-1))*(P1'*P2*P3);

e = tildeYn - Gn*hat_theta;
hat_beta = hat_theta(end);

%Auxiliary variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tildeZn = Mn(:,1:end-1);
tildeXn = Mn(:,end);

A = tildeZn*(tildeZn'*tildeZn)^(-1);
Xdd = tildeXn - A*(tildeZn'*tildeXn);
Ydd = tildeYn - A*(tildeZn'*tildeYn);
Gdd = Gn(:,end) - A*(tildeZn'*Gn(:,end));
Xddd = ((ln'*ln)^(-1))*(ln'*Xdd);

%First stage coefficient
hat_thetaFS = ((Mn'*Mn)^(-1))*(Mn'*Gdd);
hatpi = hat_thetaFS(end);
%e_AKM = Ydd - hat_beta*Gdd;

%% Compute SEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch AKMtype
    case 1
        %% AKM
        if isempty(sec_cluster_vec) == 1
            R = ( e'*ln ).^2;
            LambdaAKM = R*( Xddd.^2 );
        else
            cluster_g = unique(sec_cluster_vec);
            Nc = length(cluster_g);
            LambdaAKM_cluster = zeros(Nc,1);

            for c = 1:Nc
                lnCluster = ln(:,sec_cluster_vec'==cluster_g(c));
                XdddCluster = Xddd(sec_cluster_vec==cluster_g(c));

                RXCluster = (( e'*lnCluster  )').*XdddCluster;
                LambdaAKM_cluster(c) = sum(sum(RXCluster*RXCluster'));
            end
            LambdaAKM = sum(LambdaAKM_cluster);
        end

        %Variance matrix
        var = (1/hatpi(end)^2)*((Xdd'*Xdd)^(-1))*LambdaAKM*((Xdd'*Xdd)^(-1));
        SE = diag(var.^(1/2));

        %Confidence Interval
        CIl = hat_beta - critical*SE;
        CIu = hat_beta + critical*SE;
        CItype = 0;

        %pvalue
        tstat = (hat_beta-beta0)/SE;
        pvalue = 2*(1 - normcdf(abs(tstat),0,1));

        disp('    Coef.     | SE    | pvalue    | CIl   | CIu | CI type')
        disp([hat_beta, SE, pvalue, CIl, CIu, CItype])
    case 0
        %% AKM0
        e_null = Ydd - Gdd*beta0;

        if isempty(sec_cluster_vec) == 1
            R = ( e_null'*ln ).^2;
            LambdaAKM = R*( Xddd.^2 );
        else
            cluster_g = unique(sec_cluster_vec);
            Nc = length(cluster_g);
            LambdaAKM_cluster = zeros(Nc,1);

            for c = 1:Nc
                lnCluster = ln(:,sec_cluster_vec'==cluster_g(c));
                XdddCluster = Xddd(sec_cluster_vec==cluster_g(c));

                RXCluster = (( e_null'*lnCluster  )').*XdddCluster;
                LambdaAKM_cluster(c) = sum(sum(RXCluster*RXCluster'));
            end
            LambdaAKM = sum(LambdaAKM_cluster);
        end

        %Variance matrix
        var = (1/hatpi(end)^2)*((Xdd'*Xdd)^(-1))*LambdaAKM*((Xdd'*Xdd)^(-1));
        SE_AKMnull_n = diag(var.^(1/2));

        %pvalue
        tstat = (hat_beta-beta0)/SE_AKMnull_n;
        pvalue = 2*(1 - normcdf(abs(tstat),0,1));

        %Confidence Interval
        critical2 = critical^2;
        RY = Xdd'*Ydd;
        RX = Xdd'*Gdd;

        if isempty(sec_cluster_vec) == 1

            lnY =  Ydd'*ln ;
            lnX =  Gdd'*ln ;

            SXY = (lnY.*lnX)*( Xddd.^2 );
            SXX = (lnX.*lnX)*( Xddd.^2 );
            SYY = (lnY.*lnY)*( Xddd.^2 );

        else
            cluster_g = unique(sec_cluster_vec);
            Nc = length(cluster_g);
            SXY_cluster = zeros(Nc,1);
            SXX_cluster = zeros(Nc,1);
            SYY_cluster = zeros(Nc,1);

            for c = 1:Nc
                lnCluster = ln(:,sec_cluster_vec'==cluster_g(c));
                XdddCluster = Xddd(sec_cluster_vec==cluster_g(c));

                lnYCluster = (( Ydd'*lnCluster  )').*XdddCluster;
                lnXCluster = (( Gdd'*lnCluster  )').*XdddCluster;

                SXY_cluster(c) = sum(sum(lnYCluster*lnXCluster'));
                SXX_cluster(c) = sum(sum(lnXCluster*lnXCluster'));
                SYY_cluster(c) = sum(sum(lnYCluster*lnYCluster'));
            end
            SXY = sum(SXY_cluster);
            SXX = sum(SXX_cluster);
            SYY = sum(SYY_cluster);

        end

        Q = (RX^2)/critical2 - SXX;
        Delta = (RY*RX - critical2*SXY)^2 - (RX^2 - critical2*SXX)*(RY^2 - critical2*SYY);

        if Q > 0
            CIl = ( (RY*RX - critical2*SXY) - Delta^(1/2) )/(RX^2 - critical2*SXX);
            CIu = ( (RY*RX - critical2*SXY) + Delta^(1/2) )/(RX^2 - critical2*SXX);
            CItype = 1;
        else
            if Delta > 0
                CIl = ( (RY*RX - critical2*SXY) + Delta^(1/2) )/(RX^2 - critical2*SXX);
                CIu = ( (RY*RX - critical2*SXY) - Delta^(1/2) )/(RX^2 - critical2*SXX);
                CItype = 2;
            else
                CIl = - Inf;
                CIu =  Inf;
                CItype = 3;
            end
        end

        SE = (CIu - CIl)/(2*critical);

        disp('    Coef.     | SE    | pvalue    | CIl   | CIu | CI type')
        disp([hat_beta, SE, pvalue, CIl, CIu, CItype])


    otherwise
    disp('Specify CI type')
end

end
