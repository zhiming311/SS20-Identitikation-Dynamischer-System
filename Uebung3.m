%%
% Modellierung und Identifikation dynamischer Systeme
% Dritte Uebung: Vergleich dreier Reduktionsverfahren
% Name: Zhiming Ma
% Matrikelnummer: 3495421
% Email: zhiming0405@hotmail.com
%% Aufgabe 3.1
% x_dot = Ax + Bu
% y = Cx
% A = [-9 19   0  2 ;s
%      0  -1/2 0  0 ;
%      8  0    -1 16;
%      0 -19/2 0 -10];
% B = [1/2 5/2;
%      1/4 -1/4;
%      -1  -2;
%      3/4 1/4];
% C = [3/2 2 1 6;
%     2 3/2 -1 5];
load('ueb3_Reaktionsnetzwerk.mat');
[V,D] = eig(A);
A_hat = D;
B_hat = V^(-1)*B;
C_hat = C*V;
lambdas = diag(D);
num = length(lambdas);
ins = length(B(1,:));
outs = length(C(:,1));
disp('Diagonalisierte System ist')
disp('A_hat');
disp(A_hat);
disp('B_hat');
disp(B_hat);
disp('C_hat');
disp(C_hat);

%% stationaere Werte u_i0 = 1
g_uned = -C_hat * D^(-1) * B_hat;
y = g_uned;
mu = max(abs(g_uned),[],2);
disp('Stationaerer Ausgang y_stat ist')
disp(y);

%% Dominanzmasse
C_normiert = C_hat./mu;
B_normiert = B_hat;
D_norm = abs(C_normiert * D^(-1) * B_normiert);

D_jik_norm = zeros(outs, ins, num);
S = zeros(num,1); %Summen-Dominanzmass
M = zeros(num,1); %Maximalmass
for i = 1:num
    D_jik_norm(:,:,i) = abs(C_normiert(:,i) * B_normiert(i,:)./lambdas(i));
    S(i) = sum(D_jik_norm(:,:,i),'all');
    M(i) = max(D_jik_norm(:,:,i),[],'all');
end
disp('Summenmass S_k ist')
disp(S);
disp('Maximalmass M_k ist')
disp(M);

%% Umsotieren
um = 3; % the number of last smallest values
[minM, idx] = mink(M,um);
minEigen = lambdas(idx);
A_2_hat = diag(lambdas(idx));
A_1_hat = D;
A_1_hat(idx,:) = [];
A_1_hat(:,idx) = [];


B_2_hat = B_hat(idx,:);
B_1_hat = B_hat;
B_1_hat(idx,:) = [];

C_2_hat = C_hat(:,idx);
C_1_hat = C_hat;
C_1_hat(:,idx) = [];

num_new = length(A_1_hat(:,1));
lambdas_2 = lambdas(idx);
lambdas_1 = lambdas;
lambdas_1(idx) = [];
disp('Sortiert System ist')
disp('A_1_hat ist')
disp(A_1_hat);
disp('A_2_hat ist')
disp(A_2_hat);
disp('B_1_hat ist')
disp(B_1_hat);
disp('B_2_hat ist')
disp(B_2_hat);

%% Rekonstruktionsmatrix
Q_u =diag(ones(ins,1));
lambdas_matrix = repmat(lambdas_1, 1, num_new);
lambdas_matrix_1 = lambdas_matrix + lambdas_matrix';
lambdas_matrix_2 = repmat(lambdas_2, 1, num_new) + lambdas_matrix';

B_11 = -(B_1_hat * Q_u^2 * B_1_hat')./lambdas_matrix_1;
B_21 = -(B_2_hat * Q_u^2 * B_1_hat')./lambdas_matrix_2';

E = A_2_hat^(-1)*(B_21 + (B_2_hat - B_21*B_11^(-1)*B_1_hat)*(B_1_hat' * B_11^(-1)*B_1_hat)^(-1)*B_1_hat')*B_11^(-1)*A_1_hat;
C_reduziert = C_1_hat + C_2_hat * E;
disp('Reduzierte System :')
disp('Rekonstruktionsmatrix E ist')
disp(E);
disp('Reduzierter Ausgang C_red ist')
disp(C_reduziert);

%% Verfahren nach Guth
A_hat_strt = blkdiag(A_1_hat,A_2_hat);
B_hat_strt = [B_1_hat;B_2_hat];
C_hat_strt = [C_1_hat,C_2_hat];
L = -A_hat_strt^(-1)*B_hat_strt;
L_1p = pinv(L(1:num_new,:));
L_2 = L(num_new+1:end,:);
E_guth =L_2 * L_1p;
C_hat_guth = C_1_hat+C_2_hat*E_guth;

disp('Verfahren nach Guth')
disp('E_Guth ist')
disp(E_guth);
disp('C_1_hat_guth ist')
disp(C_hat_guth);

%% Sprungsantwort und Bodediagramm
sys = ss(A_hat_strt, B_hat_strt, C_hat_strt,0);
step(sys);
bode(sys);

%% Aufgabe 3.2
%% Ordnungsreduktion mittels balancierter Darstellung
P = lyap(A,B*B');
Q = lyap(A',C'*C);
disp('P ist')
disp(P)
disp('Q ist');
disp(Q)

%% Cholesky Zerlegung
R = chol(Q);

%% Matrix der Rechts-Eigenvekoren V_rec
[V_rec,~] = eig(R*P*R');
Sigma_sq = (V_rec'*R*P*R'*V_rec);
Sigma = diag(diag(Sigma_sq)).^(0.5);

%% transformationsmatrix
T = Sigma^(-0.5)*V'*R;
A_hat_bal = T*A*T^(-1);
B_hat_bal = T*B;
C_hat_bal = C*T^(-1);
% bis jetzt nicht richtig

%% Aufgabe 3.3