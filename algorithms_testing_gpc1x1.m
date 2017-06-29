clearvars
close all

%% Algorytm GPC 1x1 (benchmark)
% Obiekt regulacji
Gs = tf(1,[.7 1]);
Gz = c2d(Gs,0.005);
A = Gz.Denominator{1}(2);
B = Gz.Numerator{1}(2);
a = Gz.Denominator{1}(2:end);
b = Gz.Numerator{1}(2:end);
na= length(a);
nb= length(b);

% Ograniczenia
umax =  1;
umin = -1;

% Horyzonty predykcji i sterowania
N = 50; 
Nu = 50;

% Wartoœci trajektorii zadanej
yzad(   1:2000) =  .0;
yzad( 400: end) = -.1;
yzad( 800: end) =  .3;
yzad(1200: end) =  .0;
yzad(1800: end) = -.3;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1;
kk = length(yzad);

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu);
Psi = eye(N);

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(kk,1);
y = zeros(kk,1);

%% Symulacja
for k = kp:kk
    y(k) = -A*y(k-1)+B*u(k-1); % symulacja obiektu regulacji
    
    Yzad = ones(N,1)*yzad(k); % Yzad sta³e na horyzoncie predykcji
%     Yzad = yzad(min(k+(1:N),kk)); % Yzad zmienne na horyzoncie predykcji
    
    U = u((k-nb):(k-1)); % last nb values (from u(k-nb) to u(k-1))
    Y = y((k-na):k);     % last na+1 values (from y(k-na) to y(k))
    du = gpc_1x1(a,b,N,Nu,Lambda,Psi,U,Y,Yzad); % algorytm GPC 1x1
    
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