clear all; 
load("exp2_data.mat");  %load the partial noisy data for example 2
nl = data.nl; %noisy level 
Xd = data.Xd; Zd = data.Zd; Ud = data.Ud; % pre-collected data

%% Calculate the upper bound using diffent sample length
n = size(Xd,2);
bu = cell(1,n); au = cell(1,n); T = cell(1,n);
for i=1:n
   L = size(Xd{i},2);
   %upper bound claculation
   ls = 20;
   for j=ls:L
       Did = [Ud{i}(:,1:j);Zd{i}(:,1:j)];
       Did_inv = pinv(Did);
       Dib = Did_inv(:,1);  Dia = Did_inv(:,2:end);
       bu{i}(j-ls+1) = norm(Xd{i}(:,1:j)*Dib,2) + sqrt(j)*nl*norm(Dib,2); % Eq. 6
       au{i}(j-ls+1) = norm(Xd{i}(:,1:j)*Dia,2) + sqrt(j)*nl*norm(Dia,2); % Eq. 7
       T{i}(j-ls+1) = j; %data length
   end
end

figure(1)
subplot(2,2,1)
plot(T{1},bu{1},'-b')
xlabel('$L_1$',Interpreter='latex')
ylabel('$b_{1u}$',Interpreter='latex')
subplot(2,2,2)
plot(T{2},bu{2},'-b')
xlabel('$L_2$',Interpreter='latex')
ylabel('$b_{2u}$',Interpreter='latex')
subplot(2,2,3)
plot(T{1},au{1},'-b')
xlabel('$L_1$',Interpreter='latex')
ylabel('$a_{1u}$',Interpreter='latex')
subplot(2,2,4)
plot(T{2},au{2},'-b')
xlabel('$L_2$',Interpreter='latex')
ylabel('$a_{2u}$',Interpreter='latex')
