%%
% Modellierung und Identifikation dynamischer Systeme
% Zweite Uebung: Modale Ordnungsreduktion
% Name: Zhiming Ma
% Matrikelnummer: 3495421
% Email: zhiming0405@hotmail.com

%% Aufgabe 2.1
% x_dot = Ax + Bu
% y = Cx
A = [-9 19   0  2 ;s
     0  -1/2 0  0 ;
     8  0    -1 16;
     0 -19/2 0 -10];
B = [1/2 5/2;
     1/4 -1/4;
     -1  -2;
     3/4 1/4];
C = [3/2 2 1 6;
    2 3/2 -1 5];

[V,D] = eig(A);
A_hat = V^(-1)*A*V;
B_hat = V^(-1)*B;
C_hat = C*V;
lambdas = diag(D);
num = length(lambdas);
ins = length(B(1,:));
outs = length(C(:,1));


%% stationaere Werte u_i0 = 1
g_unend = -C_hat * D^(-1) * B_hat;
y = g_uned;
mu = max(abs(g_uned),[],2);
    
%Dominanzmasse
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

% Umsotieren
[minM, idx] = min(M);
minEigen = lambdas(idx);
A_1_hat = A_hat;
A_1_hat(idx,:) = [];
A_1_hat(:,idx) = [];
A_2_hat = diag(lambdas(idx));

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

% Rekonstruktionsmatrix
Q_u =diag(ones(ins,1));
lambdas_matrix = repmat(lambdas_1, 1, num_new);
lambdas_matrix_1 = lambdas_matrix + lambdas_matrix';
lambdas_matrix_2 = lambdas_1 + repmat(lambdas_2, num_new, 1);

B_11 = -(B_1_hat * Q_u^2 * B_1_hat')./lambdas_matrix_1;
B_21 = -(B_2_hat * Q_u^2 * B_1_hat')./lambdas_matrix_2';

E = A_2_hat^(-1)*(B_21 + (B_2_hat - B_21*B_11^(-1)*B_1_hat)*(B_1_hat' * B_11^(-1)*B_1_hat)^(-1)*B_1_hat')*B_11^(-1)*A_1_hat;
C_reduziert = C_1_hat + C_2_hat * E;

% Verfahren nach Guth
A_hat_strt = blkdiag(A_1_hat,A_2_hat);
B_hat_strt = [B_1_hat;B_2_hat];
L = -A_hat_strt^(-1)*B_hat_strt;
L_1p = pinv(L(1:num_new,:));
L_2 = L(end,:);
E_guth =L_2 * L_1p;
C_hat_guth = C_1_hat+C_2_hat*E_guth;