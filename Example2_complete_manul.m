clear all; 
load("exp2_cdata.mat");  %load the complete noisy data for example 2
nl = data.nl; %noisy level 
Xd = data.Xd; Zd = data.Zd; Ud = data.Ud; % pre-collected data
b = [3,3]';
%% Solve the MIPs
n = size(Xd,2);
for i=1:n
   L = size(Xd{i},2);
   %upper bound claculation
   Did = [Ud{i};Zd{i}];
   Did_inv = pinv(Did);
   Dib = Did_inv(:,1);  Dia = Did_inv(:,2:end);
   bu(i) = norm(Xd{i}*Dib,2) + sqrt(L)*nl*norm(Dib,2); % Eq. 6
   au(i) = norm(Xd{i}*Dia,2) + sqrt(L)*nl*norm(Dia,2); % Eq. 7
   % Directly set  $b_iU_{id}G_{im}=0.5$
   UGm(i) = 0.5/b(i); 
end

%% Calculate the control gain
c1 = 1; c2 = 1; %The constants that satisfy Assumption 2 
lam1 = 10; lam2 = 10; % $\lambda_1$ and $\lambda_2$

rho1 = c1*au(1); rho2 = c2*au(2); % Eq.8
k1 = round(bu(1)^2/2/lam1 + 1 + rho1,1) + 0.1; % Inq. 13, here we use round function to ensure the selection satisfy strict inequality
g1 = UGm(1)*k1; % Eq.12
r21 = 1+abs(g1);% Inq. 29
r22 = r21*rho2; % Inq. 30
r23 = max(abs(g1)*bu(1),g1^2*bu(1)+abs(g1)*rho1); % Inq. 31
r2 = r22 + r23; p2 = r2 + r2^2/4; % Inq. 32
k21 = bu(2)^2/2/lam2; k22 = lam1/2;  k2 = p2 + k21  + k22; % Inq. 14
k2 = round(k2,1) + 0.1; % Ensure the strict inequality holds
g2 = UGm(2)*k2;  % Inq. 12


