%% Algorytm DMC 1x1 (benchmark)
clearvars;
close all

% Obiekt regulacji
Gs = tf(1,[.7 1]);
Gz = c2d(Gs,0.005);
a = Gz.Denominator{1}(2);
b = Gz.Numerator{1}(2);
na= length(a);
nb= length(b);

% Ograniczenia
umax =  1;
umin = -1;

% OdpowiedŸ skokowa
S = step(Gz);
S = S(2:end); % usuwanie pierwszego elementu odpowiedzi skokowej -- KONIECZNE!
D = length(S);

%% Ogólne parametry algorytmu
% Horyzonty predykcji i sterowania
N = 50; 
Nu = 50;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1+1;
kk = 2000;
dk = 200;

% Wartoœci trajektorii zadanej
yzad = zeros(1,kk);
for k=dk:dk:kk
    yzad(k:end) = rand()-.5;
end

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu)*1.0;
Psi    = eye(N )*1.0;

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(1,kk);
y = zeros(1,kk);

%% Macierze wyznaczane offline
% Macierz M
M = zeros(N,Nu);
for row = 1:N
   for col = 1:Nu
        if(row-col+1 >= 1)
            M(row,col) = S(row-col+1);
        end
   end
end

% Macierz Mp
Mp = zeros(N,D-1);
for row = 1:N
   for col = 1:(D-1)
        Mp(row,col) = S(min(row+col,D)) - S(col);
   end
end

% Macierz K
K = (M'*Psi*M+Lambda)^(-1)*M';

%% Macierze dla wersji minimalistycznej algorytmu
Ke = sum(K(1,:));
Ku = K(1,:)*Mp;    
    
%% Symulacja
for k = kp:kk
    % symulacja obiektu regulacji
    y(k) = -a*y(k-1)+b*u(k-1); 
    
    % wprowadzanie zak³óceñ
    % y(k) = y(k) + (rand()-.5)/500;

    % wyznaczanie Yzad (sta³e na horyzoncie predykcji)    
    Yzad = ones(N,1)*yzad(k);
       
    % wyznaczanie Y
    Y    = ones(N,1)*y(k);
    
    % wyznaczanie dUp
    dUp = zeros(D-1,1);
    for p = 1:(D-1)
        if(k-p > 0); dUp(p) = u(k-p); end
        if(k-p-1 > 0); dUp(p) = dUp(p)-u(k-p-1); end
    end
    
    %% wyznaczenie du (klasycznie)
    Y0 = Y+Mp*dUp;
    dU = K*(Yzad-Y0);
    du_no = dU(1);
    
    %% wyznaczenie du (optymalnie)
    du = Ke*(yzad(k)-y(k)) - Ku*dUp;
    
    du_diff(k) = du - du_no;
    
    u(k) = u(k-1)+du;
    
    if(u(k)>umax); u(k) = umax; end
    if(u(k)<umin); u(k) = umin; end
end

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
figure;
plot(y'); hold on;
stairs(yzad,'k--'); hold off;
title('Wartoœci wyjœciowe i zadane w czasie');

figure;
stairs(u);
title('Wartoœci sterowania w czasie');

figure;
stairs(du_diff);
title('Wartoœci b³êdu w czasie');