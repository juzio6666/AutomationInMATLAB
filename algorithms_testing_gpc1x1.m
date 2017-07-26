%% Algorytm GPC 1x1 (benchmark)
clear all
global N a b na nb u y kp kk
close all

mu   =-0.165315345000003*0;%,-0.065108360000001]*0;
sigma= 0.001148139448443*0;%, 0.001082560427385]*0;

obiekt_losowy = 0;


%% Obiekt regulacji
if(obiekt_losowy == 0)
    inercje = 1;
    pobj = .7;
    ppobj = [pobj 1];
    for i=2:inercje
        ppobj = conv([pobj 1], ppobj); %(pobj(m,n)*s+1)^n
    end
    Gs = tf(1,ppobj);
    Tp = 0.005;
    Gz = c2d(Gs,Tp,'zoh');
    A = Gz.Denominator{1}(2);
    B = Gz.Numerator{1}(2);
    a = Gz.Denominator{1}(2:end);
    b = Gz.Numerator{1}(1:end);
    na= length(a);
    nb= length(b);
elseif(obiekt_losowy == 1)
    na = 5;
    nb = 10; 
    a = rand(1,na);
    b = rand(1,nb);
else
    load('parametry');
    len=80;
    a = -B_app(1:len);
    b = B_app((len+1):end);
    na= length(a);
    nb= length(b);
end

% Ograniczenia
umax =  Inf;%1;
umin = -Inf;%1;

%% Ogólne parametry algorytmu
% Horyzonty predykcji i sterowania
N  = 100; 
Nu = 100;

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
% OdpowiedŸ skokowa
S = zeros(1,N);
for j=1:N
    S(j) = 0;
    for i=1:min(j,nb); S(j) = S(j) + b(i); end
    for i=1:min(j-1,na); S(j) = S(j) - a(i)*S(j-i); end
end

% Macierz M
M = zeros(N,Nu);
for row = 1:N
   for col = 1:Nu
        if(row-col+1 >= 1)
            M(row,col) = S(row-col+1);
        end
   end
end

% Macierz K
K = (M'*Psi*M+Lambda)^(-1)*M'*Psi;

%% Macierze dla wersji minimalistycznej algorytmu
% Kyzad=zeros(1,1); % zakomentowane bo MATLAB marudzi
Ku=zeros(1,nb);     % j -> nb 
Ky=zeros(1,na+1);   % j -> (na+1)
% j -- dynamika sygna³u wejœciowego/wyjœciowego

% Kolejnoœæ nie jest przypadkowa!
fun_f(1,1); % inicjalizacja parametrów f
fun_g(1,1); % inicjalizacja parametrów g
fun_e(1,1); % inicjalizacja parametrów e

for i=1:nb
    for p=1:N
        Ku(i) = Ku(i)+K(1,p)*fun_e(p,i);
    end
end
for i=0:na
    for p=1:N
        Ky(i+1) = Ky(i+1)+K(1,p)*fun_f(p,i);
    end
end
Kyzad = sum(K(1,:));

%% Generacja macierzy
gpc1x1_matlab_to_C

%% Symulacja
du_diff = zeros(1,kk);
for k = kp:kk
    % symulacja obiektu regulacji
    ys(k) = -a*ys(k-(1:length(a)))'+b*u(k-(1:length(b)))'; 
    
    % wprowadzanie zak³óceñ
    % y(k) = y(k) + (rand()-.5)/500;
    y(k) = ys(k) + normrnd(mu,sigma);
    %y(k) = ys(k);
            
    % wyznaczanie wyjœcia modelu
    ym = 0;
    for i=1:nb; ym = ym + b(i)*u(k-i); end % u(k-i) -k=nb+1-> u(nb+1-i)
    for i=1:na; ym = ym - a(i)*y(k-i); end % y(k-i) -k=na+1-> y(na+1-i)
    
    % wyznaczanie d
    d = y(k) - ym;
    
    % wyznaczanie Y0
    Y0 = zeros(N,1);
    for p=1:N
        Y0(p) = d;
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
    
    % wyznaczanie Yzad (sta³e na horyzoncie predykcji)
    Yzad = ones(N,1)*yzad(k);
    
    %% wyznaczenie du (pó³optymalnie)
    du_po = K(1,:)*(Yzad-Y0);
    
    %% wyznaczenie du (optymalnie)
    du = Kyzad*yzad(k) - Ku*u(k-(1:nb))' - Ky*y(k-(0:na))';
    
    du_diff(k) = du - du_po;
    
    u(k) = u(k-1)+du;    
    
    if(u(k)>umax); u(k) = umax; end
    if(u(k)<umin); u(k) = umin; end
end

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
figure;
plot((0:(kk-1))*Tp,y'); hold on;
stairs((0:(kk-1))*Tp,yzad','k--'); hold off;
title('Wartoœci wyjœciowe i zadane w czasie');

figure;
stairs((0:(kk-1))*Tp,u');
title('Wartoœci sterowania w czasie');

% figure;
% stairs(du_diff);
% title('Wartoœci b³êdu w czasie');
csvwrite('../LaTeX/DUNNO_Pomiar_czasu_algorytmow_regulacji/dane/gpc1x1simdist.csv',[((1:2000)'-1)*Tp y' yzad' u']);

%% Funkcje do wyznaczania minimalnej postaci algorytmu GPC
function out = fun_g(p,j)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb 
    global G
    
    o = N;
    if(isempty(G))
        G=cell(N,nb-1+o); 
        for p=1:N; for j=(1-p):(nb-1); fun_g(p,j); end; end
    end
    if(~isempty(G{p,j+o}))
        out = G{p,j+o};
    else
        if(p == 1)
            if(j>=0 && j<=(nb-1))
                out = b(j+1);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(G);disp(o);error('Error!');
            end
        elseif(p>=2 && p<=N)
            if(j>=(1-p) && j<=(-1))
                out = fun_g(p-1,j+1);
            elseif(j>=0 && j<=(nb-2))
                out = fun_f(p-1,0)*fun_g(1,j)+fun_g(p-1,j+1);
            elseif(j==(nb-1))
                out = fun_f(p-1,0)*fun_g(1,j);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(G);disp(o);error('Error!');
            end
        else
            disp(N);disp(a);disp(b);disp(na);disp(nb);disp(G);disp(o);error('Error!');
        end
        G{p,j+o} = out;
    end
end


function out = fun_e(p,j)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb 
    global E
    
    o = 0;
    if(isempty(E))
        E=cell(N,nb+o); 
        for p=1:N; for j=(1-o):nb; fun_e(p,j); end; end
    end
    if(~isempty(E{p,j+o}))
        out = E{p,j+o};
    else
        if(j==1 && nb==1)
            out = 0;
        elseif(j==1 && nb>1)
            out = fun_g(p,j);
        elseif(j>=2 && j<=(nb-1) && j<nb && nb>1)
            out = fun_g(p,j) - fun_g(p,j-1);
        elseif(j==nb && nb>1)
            out = -fun_g(p,j-1);
        else
            disp(N);disp(a);disp(b);disp(na);disp(nb);disp(E);disp(o);error('Error!');
        end
        E{p,j+o} = out;
    end
end

function out = fun_f(p,j)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb 
    global F
    
    o = 1;
    if(isempty(F))
        F=cell(N,na+o);
        for p=1:N; for j=(1-o):na; fun_f(p,j); end; end
    end
    if(~isempty(F{p,j+o}))
        out = F{p,j+o};
    else
        if(p == 1)
            if(j==0)
                out = 1 - a(1);
            elseif(j>=1 && j<=(na-1))
                out = a(j)-a(j+1);
            elseif(j==na)
                out = a(j);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(F);disp(o);error('Error!');
            end
        elseif(p>=2 && p<=N)
            if(j>=0 && j<=(na-1))
                out = fun_f(p-1,0)*fun_f(1,j)+fun_f(p-1,j+1);
            elseif(j==na)
                out = fun_f(p-1,0)*fun_f(1,j);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(F);disp(o);error('Error!');
            end
        else
            disp(N);disp(a);disp(b);disp(na);disp(nb);disp(F);disp(o);error('Error!');
        end
        F{p,j+o} = out;
    end
end