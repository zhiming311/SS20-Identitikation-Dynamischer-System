%%
% Modellierung und Identifikation dynamischer Systeme
% Vierte Uebung: System Identification Toolbox
% Name: Zhiming Ma
% Matrikelnummer: 3495421
% Email: zhiming0405@hotmail.com

clear;
%% Aufgabe 4.1 Vorbereitung

%% Aufgabe 4.2 Black-Box-Identifikation
load('Uebung_4_Data\Ueb4.mat');
Chirp_iddata = iddata(dataChirp.y, dataChirp.u, dataChirp.Ts);
Valid_iddata = iddata(dataValid.y, dataValid.u, dataValid.Ts);
fChirp = figure('Name','Chirp');
plot(Chirp_iddata);
fValid = figure('Name','Valid');
plot(Valid_iddata);

% Befehl: ident (System Identification Toolbox ist benoetigt)
ident;
% Import data-> data object
% Waehlen Chirp_iddata als "working data", Valid_iddata als "Validation Data"
% Den Frequenzgang der Identifikationsdaten anzeigen

% Polynomial Structure ARX421, ARX331 sind die Besten
% Mit detrend, besser Ergebnisse sind erreichbar

%% Aufgabe 4.3 Identifikation mittels Regression und Instrumental Variables
%% a) Durch Euler-vorwaerts Verfahren
% G(z) = (5e-5)/(z^2-1.999*z+0.99915)
m = 2;
d = 0.2;
k = 3;
T_s = 0.01;

%% b) Differenzengleichung 
% y(k+2) + a_1 * y(k+1) + a_0 * y(k) = b_0 * u(k)
% G(z) = b_0 / (z^2 + a_1 * z + a_0)
%
%% c) Identifikation mittels Regression
load('Uebung_4_data\ES_iddata_no_noise.mat');
num_a = 2;
num_b = 1;
num_measure = length(y);

psi = zeros(num_measure - num_a, num_a + num_b); % num_measure - num_a sind taetig Messungen

for i = 1: num_measure - num_a
    psi(i,:) = [-y(i+num_a-1:-1:i)', u(i)'];
end

y_taetig = y(num_a+1:num_measure);
% y_schatz = psi * theta;
% err = y_taetig - y_schatz;
% J = err'*err;
theta = (psi'*psi)\psi'*y_taetig;
% theta = [-1.999, 0.9992, 0.00005]; Fast gleich wie der Tatsach

%% d) Validierung
load('Uebung_4_data\ES_validdata_no_noise.mat')
y_pre = zeros(num_measure,1);
for i = 1:num_measure - num_a
    y_pre(i+num_a) = [-y_pre(i+num_a-1), -y_pre(i+num_a-2), u(i+num_a-2)]*theta;
end

fig = figure('Name','Vergleichen zwischen simuliert Antwort und Validierungsdaten');
plot(y_pre,'r');
hold on;
plot(y,'b');
legend;

%% e) Identifikation mittels Regression mit verrauschten Messdaten
clear;
close;

load('Uebung_4_data\ES_iddata_noise.mat');
num_a = 2;
num_b = 1;
num_measure = length(y);

psi = zeros(num_measure - num_a, num_a + num_b); % num_measure - num_a sind taetig Messungen

for i = 1: num_measure - num_a
    psi(i,:) = [-y(i+num_a-1:-1:i)', u(i)'];
end

y_taetig = y(num_a+1:num_measure);
theta = (psi'*psi)\psi'*y_taetig;

y_pre = zeros(num_measure,1);
for i = 1:num_measure - num_a
    y_pre(i+num_a) = [-y_pre(i+num_a-1), -y_pre(i+num_a-2), u(i+num_a-2)]*theta;
end

load('Uebung_4_data\ES_validdata_noise.mat')
fig = figure('Name','Vergleichen zwischen simuliert Antwort und Validierungsdaten');
plot(y_pre,'r');
hold on;
plot(y,'b');
legend;


