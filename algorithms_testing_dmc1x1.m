%% Algorytm DMC 1x1 (benchmark)
clear all

close all

mu   =-0.165315345000003*0;%,-0.065108360000001]*0;
sigma= 0.001148139448443*0;%, 0.001082560427385]*0;

obiekt_losowy = 1;


%% Obiekt regulacji
Tp = 0.005;
if(obiekt_losowy == 0)
    inercje = 1;
    pobj = .7;
    ppobj = [pobj 1];
    for i=2:inercje
        ppobj = conv([pobj 1], ppobj); %(pobj(m,n)*s+1)^n
    end
    Gs = tf(1,ppobj);
    Gz = c2d(Gs,Tp,'zoh');
    A = Gz.Denominator{1}(2);
    B = Gz.Numerator{1}(2);
    a = Gz.Denominator{1}(2:end);
    b = Gz.Numerator{1}(2:end);
    na= length(a);
    nb= length(b);
elseif(obiekt_losowy == 1)
    Tp=0.01;
    K=0.9690;
    T1=0.032857;
    T2=0.331948;
    Gs=tf(K,[T1*T2,T1+T2,1]);
    Gz = c2d(Gs,Tp,'zoh');
    A = Gz.Denominator{1}(2);
    B = Gz.Numerator{1}(2);
    a = Gz.Denominator{1}(2:end);
    b = Gz.Numerator{1}(2:end);
    na= length(a);
    nb= length(b);
else
    na = 1;
    nb = 1; 
    a = rand(1,na);
    b = rand(1,nb);
    Gz = tf([0 b],[1 a], Tp);
end

% Ograniczenia
umax =  1;
umin = -1;

% OdpowiedŸ skokowa
S = step(Gz);
S = S(1:end); % usuwanie pierwszego elementu odpowiedzi skokowej -- KONIECZNE! pod warunkiem, ¿e czas regulacji jest znikomy w stosunku do okresu próbkowania
SS = S;
DD = length(S);
D = 1000;
if(DD<D)
    S(DD:D)=S(DD);
else
    S=S(1:D);
end 

%% Ogólne parametry algorytmu
% Horyzonty predykcji i sterowania
N = D; 
Nu = D;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1+1;
kk = 2000;
dk = 200;

% Wartoœci trajektorii zadanej
yzad = zeros(1,kk);
% for k=dk:dk:kk
%     yzad(k:end) = (rand()*2-1)*0.1;
% end
yzad(100:end) = -.1;
yzad(500:end) =  .1;
yzad(900:end) = -.3;
yzad(1300:end) = .2;
yzad(1700:end) = .0;

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu)*1.0;
Psi    = eye(N )*1.0;

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(1,kk);
y = zeros(1,kk);
ys = zeros(1,kk);

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
K = (M'*Psi*M+Lambda)^(-1)*M'*Psi;

%% Macierze dla wersji minimalistycznej algorytmu
Ke = sum(K(1,:));
Ku = K(1,:)*Mp;    
    
%% Generacja macierzy
dmc1x1_matlab_to_C

%% Symulacja
du_diff = zeros(1,kk);
for k = kp:kk
    ys(k) = -a*ys(k-(1:length(a)))'+b*u(k-(1:length(b)))'; 
    
    % wprowadzanie zak³óceñ
    % y(k) = y(k) + (rand()-.5)/500;
    y(k) = ys(k) + normrnd(mu,sigma);

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
stairs(u');
title('Wartoœci sterowania w czasie');

% csvwrite('../LaTeX/DUNNO_Pomiar_czasu_algorytmow_regulacji/dane/dmc1x1simdist.csv',[((1:2000)'-1)*Tp y' yzad' u']);