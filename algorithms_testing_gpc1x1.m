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

% Memmory assignment
M = zeros(N,Nu);
S = zeros(N,1);

% Vector S
for j=1:N
    S(j) = 0;
    for i=1:min(j,nb); S(j) = S(j) + b(i); end
    for i=1:min(j-1,na); S(j) = S(j) - a(i)*S(j-i); end
end

% Matrix M
for row = 1:N
   for col = 1:Nu
        if(row-col+1 >= 1)
            M(row,col) = S(row-col+1);
        end
   end
end

K = (M'*Psi*M+Lambda)^(-1)*M';    
Knu = K(1,:);

%% Symulacja
for k = kp:kk
    y(k) = -A*y(k-1)+B*u(k-1); % symulacja obiektu regulacji
    
    Yzad = ones(N,1)*yzad(k); % Yzad sta³e na horyzoncie predykcji
    
    U = u((k-nb):(k-1)); % last nb values (from u(k-nb) to u(k-1))
    Y = y((k-na):k);     % last na+1 values (from y(k-na) to y(k))
    
    
    ym = 0;
    for i=1:nb; ym = ym + b(i)*u(k-i); end % u(k-i) -k=nb+1-> u(nb+1-i)
    for i=1:na; ym = ym - a(i)*y(k-i); end % y(k-i) -k=na+1-> y(na+1-i)
    d = y(k) - ym;
    
    Y0 = ones(N,1)*d;
    for p=1:N
        for i=1:nb
            if( p-i <= -1 )
                Y0(p) = Y0(p) + b(i)*u(k+p-i); % u(k+p-i) -k=nb+1-> u(nb+1+p-i)
            else
                Y0(p) = Y0(p) + b(i)*u(k-1);   % u(k-1) -k=nb+1-> u(nb+1-1)
            end 
        end
        for i=1:na
            if( p-i <= 0 )
                Y0(p) = Y0(p) - a(i)*y(k+p-i); % y(k+p-i) -k=na+1-> y(na+1+p-i)        
            else
                Y0(p) = Y0(p) - a(i)*Y0(p-i);       % Y0(k+p-i) -k=0-> Y0(p-i)
            end
        end 
    end
    
    du = Knu*(Yzad-Y0);
    
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