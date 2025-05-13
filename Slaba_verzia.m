clc
clear
close all

format long

%diferencialne rovnice Van der Vusse CSTR
% dCA/dt=-k1*CA-k3*CA^2+F/V*(CAf-CA)
% dCB/dt=k1*CA-k2*CB+F/V*CB
% y=CB

% parametre chemickeho reaktora
k1=50;
k2=100;
k3=10;
V=1;
% CAf=6;

% F0=1;
% F0=52;

CAf=input('Zadaj CAf [mol/l]:')  %z int. <1;10>52
F0=input('Zadaj F0 [l/h]:')     %z int. <0.1;60>

%volba pracovneho bodu, vypocet ustalenych stavov CA0 a CB0 pre vstup F0
CA01=(-k1-F0/V+sqrt((k1+F0/V)^2+4*k3*F0/V*CAf))/(2*k3);
CA02=(-k1-F0/V-sqrt((k1+F0/V)^2+4*k3*F0/V*CAf))/(2*k3);

if CA01>0    
    CA0=CA01;
else
    CA0=CA02;
end

sprintf('Parameter CA0 je: %f',CA0)

CB0=k1*CA0/(k2+F0/V);

sprintf('Parameter CB0 je: %f',CB0)

%ziskanie prenosovej funkcie pre vybrany PB
a=-k1-2*k3*CA0-F0/V;
b=CAf/V-CA0/V;
l=-k2-F0/V;
m=-CB0/V;

disp('Prenosova funkcia je:')
g= tf ([m -a*m+b*k1],[1 -a-l a*l])

%value of denominator
sNum = g.Numerator{1};
sDen = g.Denominator{1};

K = ( (sNum(1, 3)) / (sDen(1, 3)) );

%% Určenie času stabilizácie a periódy vzorkovania
info = stepinfo(g);
t_ust = info.SettlingTime; % Čas nastavenia systému
T = ((1/10) * t_ust); %Perioda vzorkovania 


%% Diskretizácia povodnej prenosovej funkcie g
g_z = c2d(g, T); % g_z - diskretna prenosova charakteristika g.

B_z = g_z.Numerator; % B(z) - citatel
A_z = g_z.Denominator; % A(z) - menovatel

syms z;
integrator = poly2sym([-1, 1], z);

%% Stabilita A_z

A_z_roots = roots(A_z{1});

roots_stable = A_z_roots(abs(A_z_roots) < 1);      %stabilne
roots_unstable = A_z_roots(abs(A_z_roots) >= 1); 


A_plus = poly(roots_stable);       % A^+(z)
A_minus = poly(1 ./ roots_unstable); % A^-(z)

%Control
% A_rec = conv(A_plus, A_minus);  % Должен получиться оригинальный полином (почти)
% disp('A(z) оригинальный:'), disp(A_z);
% disp('A(z) восстановленный:'), disp(A_rec);

A = A_plus * A_minus;


%% Stabilita B_z

B_z_roots = roots(B_z{1});

roots_B_stable = B_z_roots(abs(B_z_roots) < 1);      %stabilne
roots_B_unstable = B_z_roots(abs(B_z_roots) >= 1); 

B_plus = poly(roots_B_stable);       % B^+(z) = 1 (если нет стабильных корней)
B_minus = poly(roots_B_unstable);  

% Если B^-(z) должен включать исходный полином (по аналогии с примером из скрина)
if isempty(roots_B_unstable)
    B_minus = 1;  % Если все корни стабильны (но в вашем случае это не так)
else
    B_minus = B_z{1};  % Используем исходный B(z), так как все его корни нестабильны
end

% Control
% disp('Стабильная часть B^+(z):'); disp(B_plus);
% disp('Нестабильная часть B^-(z):'); disp(B_minus);

B = B_plus * B_minus;




%% Diof rovnica

syms z p0_0 p0_1 q0_0;
B_minus_poly = vpa((B_minus(1,2)*z + B_minus(1, 3) *z^2),6);
P0_z = poly2sym([p0_1, p0_0], z);
Q0_z = q0_0;


term1 = expand(integrator * P0_z);
simplified_term1 = collect(term1, z);  %Control

term2 = expand(B_minus_poly * Q0_z);
simplified_term2 = vpa(collect(term2, z), 6); %Control

right_side = 0.0016;


left_expr = expand(integrator * P0_z + B_minus_poly * Q0_z);

%%lava strana upravena
res_expr = vpa(collect(left_expr, z), 6); %Control left side good form



%% Output coefficients 

[coeffs, terms] = coeffs(res_expr, z);
% Коэффициенты для разных степеней z
coeff_z_2 = vpa(coeffs(terms == z^2), 6);  %  z^-2
coeff_z_1 = vpa(coeffs(terms == z^1), 6);  %  z^-1
coeff_z_0 = vpa(coeffs(terms == z^0), 6); %   z^0

%% urobime rovnice
eq1 = coeff_z_2 == 0;      
eq2 = coeff_z_1 == 0;      
eq3 = coeff_z_0 == right_side; 

%% Finding coeffs p0_0, p0_1, q0_0

sol = solve([eq1, eq2, eq3], [p0_0, p0_1, q0_0]);

% coeffs result
% fprintf('Решение системы:\n');
% fprintf('p0_0 = %.8f\n', double(sol.p0_0));
% fprintf('p0_1 = %.8f\n', double(sol.p0_1));
% fprintf('q0_0 = %.8f\n', double(sol.q0_0));


% Записываем коэффициенты в переменные
p0_0_value = double(sol.p0_0); 
p0_1_value = double(sol.p0_1); 
q0_0_value = double(sol.q0_0);  


%% P(z) a Q(z)

syms z;
Q_z = (q0_0_value) * A_plus; %numerator
P_z = conv([1 -1], [p0_0_value p0_1_value] ); %denominator

%% %% Getting data for output

simOut = sim("schema_slaba.slx");


oneDOFdata = simOut.get('Slaba_verzia_1Dof'); % get data from Scope
oneDOFdata_time = oneDOFdata.time; % Data of time
oneDOFScopeData_sig1 = oneDOFdata.signals(1).values; 
oneDOFScopeData_sig2 = oneDOFdata.signals(2).values; 


%Data for 1 DOF model Akcny zasah
twoDOFAkZasdata = simOut.get('Slaba_verzia_baza'); % get data from Scope
twoDOFAkZasScopeData_time = twoDOFAkZasdata.time; % Data of time
twoDOFAkZasScopeData_sig1 = twoDOFAkZasdata.signals.values(:,1); 



%% Output results 

% 1 DOF model URO
subplot(2,2,1);
plot(oneDOFdata_time, oneDOFScopeData_sig1, 'k', 'LineWidth', 1, 'LineStyle', '--'); 
hold on;
plot(oneDOFdata_time, oneDOFScopeData_sig2(:,1), 'r', 'LineWidth', 1); 
xlabel('Čas  (s)'); 
ylabel('Amplitúda'); 
title('Slaba verzia - 1 DOF'); 
grid on; 

% 1 DOF model Akcny Zasah
subplot(2,2,2);
stairs(twoDOFAkZasScopeData_time, twoDOFAkZasScopeData_sig1, 'b', 'LineWidth', 1);
xlabel('Čas  (s)'); 
ylabel('Amplitúda'); 
title('Slaba verzia baza');
grid on;

