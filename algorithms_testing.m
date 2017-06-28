%clearvars
%close all

%% Algorytm DMC 1x1 (benchmark)
% Obiekt regulacji
Gs = tf(1,[.7 1]);
Gz = c2d(Gs,0.005);
a = Gz.Denominator{1}(2);
b = Gz.Numerator{1}(2);

umax =  1;
umin = -1;

% OdpowiedŸ skokowa
S = step(Gz);
S = S(2:end); % usuwanie pierwszego elementu odpowiedzi skokowej -- KONIECZNE!

% Horyzonty predykcji i sterowania
N = 50; 
Nu = 50;

% Wartoœci trajektorii zadanej
yzad(1:2000) = 0;
yzad(   1: end) = -.8;
%yzad( 200: end) =  .1;
yzad( 400: end) = -.1;
%yzad( 600: end) =  .7;
yzad( 800: end) = -.3;
%yzad(1000: end) = -.4;
yzad(1200: end) =  .0;
%yzad(1400: end) =  .1;
%yzad(1600: end) =  .9;
yzad(1800: end) = -.3;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = 2;
kk = length(yzad);

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu);
Psi = eye(N);

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
%u = zeros(kk,1);
%y = zeros(kk,1);

% Horyzont dynamiki (do wyznaczania dUp)
D = length(S);

%% Symulacja
for k = kp:kk
    y(k) = -a*y(k-1)+b*u(k-1); % symulacja obiektu regulacji
    
    Yzad = ones(N,1)*yzad(k); % Yzad sta³e na horyzoncie predykcji
%     Yzad = yzad(min(k+(1:N),kk)); % Yzad zmienne na horyzoncie predykcji
    Y = y(k);
    dUp = zeros(D-1,1);
    for p = 1:(D-1)
        if(k-p > 0); dUp(p) = u(k-p); end
        if(k-p-1 > 0); dUp(p) = dUp(p)-u(k-p-1); end
    end
    
    du = dmc_1x1(S,N,Nu,Lambda,Psi,dUp,Y,Yzad); % algorytm DMC 1x1
    
    u(k) = u(k-1)+du;
    if(u(k)>umax); u(k) = umax; end
    if(u(k)<umin); u(k) = umin; end
end

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
figure;
plot(y); hold on;
stairs(yzad,'k--'); hold off;

figure;
stairs(u);