%% Algorytm DMC 2x2 (benchmark)
clearvars
global N a b na nb nu ny
close all

%% Obiekt regulacji
%      input 1        input 2
Gs = [tf( 1,[.7 1]), tf( 5,[.03 1]);  % output 1
      tf( 1,[.05 1]), tf( 2,[.4 1])]; % output 2
Tp = 0.05;
Gz = c2d(Gs,Tp,'zoh');
ny = 2;
nu = 2;

% Based on STP script, pages 79 (tab. 3.1) and 86 (eq. 3.49)
A = zeros(2,2);
B = zeros(2,2);
for m=1:2
    for n=1:2
        alfa = exp(-1/Gs(m,n).Denominator{1}(1)*Tp);
        B(m,n) = Gs(m,n).Numerator{1}(2)*(1-alfa);
        A(m,n) =                           -alfa ;
    end
end

% Y1/U1=G(1,1) and Y1/U2=G(1,2) => Y1 = G(1,1)*U1 + G(1,2)*U2
for m=1:2
    a(m,:)   = [A(m,1)+A(m,2), A(m,1)*A(m,2)];
    b(m,1,:) = [0, B(m,1), A(m,2)*B(m,1)];
    b(m,2,:) = [0, B(m,2), A(m,1)*B(m,2)];
end

na = size(a,2);
nb = size(b,3);

umax =  1;
umin = -1;

%% Ogólne parametry algorytmu
% Horyzonty predykcji i sterowania
N  = 50; 
Nu = 50;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1+1;
kk = 2000;
dk = 200;

% Wartoœci trajektorii zadanej
yzad = zeros(ny,kk);
for k=dk:dk:kk
    for m=1:ny
        yzad(m,(k-(m-1)*dk/ny):end) = rand()-.5;
    end
end

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu*nu)*1.0;
Psi    = eye(N *ny)*1.0;

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(nu,kk);
y = zeros(ny,kk);

% OdpowiedŸ skokowa
S = zeros(ny,nu,N);
for k = 1:size(S,3)
    for m=1:ny
        for n=1:nu
            for i=1:min(k,nb)
                S(m,n,k) = S(m,n,k) + b(m,n,i)*1;
            end
            for i=1:min(k-1,na)
                S(m,n,k) = S(m,n,k) - a(m,i)*S(m,n,k-i);
            end         
        end   
    end 
end

% Macierz M
M = cell(N,Nu);
for row = 1:N
   for col = 1:Nu
        if(row-col+1 >= 1)
            M{row,col} = S(:,:,row-col+1);
        else
            M{row,col} = zeros(size(S(:,:,1)));
        end
   end
end
M=cell2mat(M);

%% Macierze wyznaczane offline
K = (M'*Psi*M+Lambda)^(-1)*M';

%% Macierze dla wersji minimalistycznej algorytmu
Kyzad = zeros(nu,ny);
Ku = zeros(nu,nu,nb);   % r,n,j -> nu x nu x nb
Ky = zeros(nu,ny,na+1); % r,m,j -> nu x ny x (na+1)
% r -- numer sygna³u steruj¹cego, którego przyrost jest wyliczany
% n -- numer sygna³u steruj¹cego
% m -- numer sygna³u wyjœciowego
% j -- dynamika sygna³u wejœciowego/wyjœciowego

% Kolejnoœæ nie jest przypadkowa!
fun_f(1,1,1);   % inicjalizacja parametrów f
fun_g(1,1,1,1); % inicjalizacja parametrów g
fun_e(1,1,1,1); % inicjalizacja parametrów e

%% Wyznaczanie parametrów Ku, Ky
for r=1:nu
    for n=1:nu
        for j=1:nb
            for p=1:N
                for m=1:ny
                    s=(p-1)*ny+m;
                    Ku(r,n,j) = Ku(r,n,j) - K(r,s)*fun_e(p,j,m,n);
                end
            end
        end
    end
end
for r=1:nu
    for m=1:ny
        for j=0:na
            for p=1:N
                s=(p-1)*ny+m;
                Ky(r,m,j+1) = Ky(r,m,j+1) - K(r,s)*fun_f(p,j,m);
            end
        end
    end
end
for r=1:nu
    for m=1:ny
        for p=1:N
            s=(p-1)*ny+m;
            Kyzad(r,m) = Kyzad(r,m) + K(r,s);
        end
    end
end

%% Symulacja
dudiff = zeros(nu,kk);
for k = kp:kk
    % symulacja obiektu regulacji
    for m=1:ny
        for n=1:nu
            for i=1:nb
                if(k-i>=1)
                    y(m,k) = y(m,k) + b(m,n,i)*u(n,k-i);
                end
            end
        end
        for i=1:na
            if(k-i>=1)
                y(m,k) = y(m,k) - a(m,i)*y(m,k-i);
            end
        end         
    end 
    
    % wprowadzanie zak³óceñ
    % for m=1:ny; y(m,k) = y(m,k) + (rand()-.5)/500; end;
     
    % wyznaczanie wyjœcia modelu
    ym = zeros(ny,1);
    for m=1:ny
        for n=1:nu
            for i=1:nb
                if(k-i>=1)
                    ym(m) = ym(m) + b(m,n,i)*u(n,k-i);
                end
            end
        end
        for i=1:na
            if(k-i>=1)
                ym(m) = ym(m) - a(m,i)*y(m,k-i);
            end
        end            
    end 
    
    % wyznaczanie d
    d = y(:,k)-ym;
    
    % wyznaczanie Y0
    Y0=zeros(ny,N);
    for m=1:ny
        for p=1:N
            Y0(m,p) = d(m);
            for n=1:nu
                for i=1:nb
                    if(-i+p<=-1)
                        Y0(m,p) = Y0(m,p) + b(m,n,i)*u(n,k-i+p);
                    else
                        Y0(m,p) = Y0(m,p) + b(m,n,i)*u(n,k-1);
                    end
                end
            end
            for i=1:na
                if(-i+p<=0)
                    Y0(m,p) = Y0(m,p) - a(m,i)*y(m,k-i+p);
                else
                    Y0(m,p) = Y0(m,p) - a(m,i)*Y0(m,-i+p);
                end
            end            
        end
    end
    Y0 = reshape(Y0,[],1);
    
    % wyznaczanie Yzad (sta³e na horyzoncie predykcji)
    Yzad = repmat(eye(ny),N,1)*yzad(:,k); 
    
    %% wyznaczenie du (pó³optymalnie)
    du_po = K(1:nu,:)*(Yzad-Y0);
    
    %% wyznaczenie du (optymalnie)
    du = zeros(nu,1);
    for r=1:nu
        for m=1:ny
            du(r) = du(r) + Kyzad(r,m)*yzad(m,k);
        end
        for n=1:nu
            for j=1:nb
                du(r) = du(r) + Ku(r,n,j)*u(n,k-j);
            end
        end
        for m=1:ny
            for j=0:na
                du(r) = du(r) + Ky(r,m,j+1)*y(m,k-j);
            end
        end
    end
    
    du_diff(1:nu,k) = du-du_po;
    
    u(:,k) = u(:,k-1)+du;
    
    for n=1:nu
        if(u(n,k)>umax); u(n,k) = umax; end
        if(u(n,k)<umin); u(n,k) = umin; end
    end
end

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
figure;
plot(y'); hold on;
stairs(yzad','k--'); hold off;
title('Wartoœci wyjœciowe i zadane w czasie');

figure;
stairs(u');
title('Wartoœci sterowania w czasie');

figure;
stairs(du_diff');
title('Wartoœci b³êdu w czasie');

%% Funkcje do wyznaczania minimalnej postaci algorytmu GPC
function out = fun_e(p,j,m,n)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb nu ny 
    persistent E o
    
    if(isempty(E))
        o = 0;
        E=cell(ny,nu,N,nb+o);
        for m=1:ny; for n=1:nu; for p=1:N; for j=(1-o):nb; fun_e(p,j,m,n); end; end; end; end
    end
    if(~isempty(E{m,n,p,j+o}))
        out = E{m,n,p,j+o};
    else
        if(j==1 && nb==1)
            out = 0;
        elseif(j==1 && nb>1)
            out = fun_g(p,j,m,n);
        elseif(j>=2 && j<=(nb-1) && j<nb && nb>1)
            out = fun_g(p,j,m,n) - fun_g(p,j-1,m,n);
        elseif(j==nb && nb>1)
            out = -fun_g(p,j-1,m,n);
        else
            disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);
            disp(p);disp(j);disp(m);disp(n);
            %disp(E);disp(o);
            error('Error!');
        end
        E{m,n,p,j+o} = out;
    end
end

function out = fun_f(p,j,m)
    global N a b na nb ny nu
    persistent F o
    
    if(isempty(F))
        o = 1;
        F=cell(ny,N,na+o);
        for m=1:ny; for p=1:N; for j=(1-o):na; fun_f(p,j,m); end; end; end
    end
    if(~isempty(F{m,p,j+o}))
        out = F{m,p,j+o};
    else
        if(p == 1)
            if(j==0)
                out = 1 - a(m,1);
            elseif(j>=1 && j<=(na-1))
                out = a(m,j)-a(m,j+1);
            elseif(j==na)
                out = a(m,j);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);disp(F);disp(o);error('Error!');
            end
        elseif(p>=2 && p<=N)
            if(j>=0 && j<=(na-1))
                out = fun_f(p-1,0,m)*fun_f(1,j,m)+fun_f(p-1,j+1,m);
            elseif(j==na)
                out = fun_f(p-1,0,m)*fun_f(1,j,m);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);disp(F);disp(o);error('Error!');
            end
        else
            disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);disp(F);disp(o);error('Error!');
        end
        F{m,p,j+o} = out;
    end
end

function out = fun_g(p,j,m,n)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb nu ny
    persistent G o
    
    if(isempty(G))
        o = N;
        G=cell(ny,nu,N,nb-1+o);
        for m=1:ny; for n=1:nu; for p=1:N; for j=(1-p):(nb-1); fun_g(p,j,m,n); end; end; end; end
    end
    if(~isempty(G{m,n,p,j+o}))
        out = G{m,n,p,j+o};
    else
        if(p == 1)
            if(j>=0 && j<=(nb-1))
                out = b(m,n,j+1);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);disp(G);disp(o);error('Error!');
            end
        elseif(p>=2 && p<=N)
            if(j>=(1-p) && j<=(-1))
                out = fun_g(p-1,j+1,m,n);
            elseif(j>=0 && j<=(nb-2))
                out = fun_f(p-1,0,m)*fun_g(1,j,m,n)+fun_g(p-1,j+1,m,n);
            elseif(j==(nb-1))
                out = fun_f(p-1,0,m)*fun_g(1,j,m,n);
            else
                disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);disp(G);disp(o);error('Error!');
            end
        else
            disp(N);disp(a);disp(b);disp(na);disp(nb);disp(nu);disp(ny);disp(G);disp(o);error('Error!');
        end
        G{m,n,p,j+o} = out;
    end
end