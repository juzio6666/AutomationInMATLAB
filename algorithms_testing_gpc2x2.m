%% Algorytm DMC 2x2 (benchmark)
clear all
global N a b na nb nu ny
%close all

mu   =[-0.165315345000003,-0.065108360000001]*0;
sigma=[ 0.001148139448443, 0.001082560427385]*0;

obiekt_losowy = 0;
ny = 2;
nu = 2;
%% Obiekt regulacji
if(obiekt_losowy == 0)
    inercje = 2;

    pobj = [.7, .3; .5, .4];
    ppobj = cell(2,2);
    for m=1:2
        for n=1:2
            ppobj{m,n} = [pobj(m,n) 1];
            for i=2:inercje
                ppobj{m,n} = conv([pobj(m,n) 1], ppobj{m,n}); %(pobj(m,n)*s+1)^n
            end
        end
    end

    %     input 1            input 2
    Gs = [tf( 1,ppobj{1,1}), tf( 5,ppobj{1,2});  % output 1
          tf( 1,ppobj{2,1}), tf( 2,ppobj{2,2})]; % output 2
    Tp = 0.005;
    Gz = c2d(Gs,Tp,'zoh');

    % Y1/U1=G(1,1) and Y1/U2=G(1,2) => Y1 = G(1,1)*U1 + G(1,2)*U2
    for m=1:2
        tmpa = conv(Gz.Denominator{m,1},Gz.Denominator{m,2});
        a(m,:) = tmpa(2:end);
    
        tmpb = conv(Gz.Numerator{m,1},Gz.Denominator{m,2});
        b(m,1,:) = tmpb;
    
        tmpb = conv(Gz.Numerator{m,2},Gz.Denominator{m,1});
        b(m,2,:) = tmpb;
    end
    na = size(a,2);
    nb = size(b,3);
else
    na = 10;
    nb = 30; 
    for m=1:ny
        a(m,:)   = rand(1,na);
        for n=1:nu
            if(m==n)
                b(m,n,:) = rand(1,1,nb);
            else
                b(m,n,:) = rand(1,1,nb)*0.1;
            end
        end
    end
end

% Ograniczenia
umax =  1;
umin = -1;

%% Ogólne parametry algorytmu
% Horyzonty predykcji i sterowania
N  = 500; 
Nu = 50;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1+1;
kk = 2000;
dk = 200;

% Wartoœci trajektorii zadanej
yzad = zeros(ny,kk);
for k=dk:dk:kk
    for m=1:ny
        yzad(m,(k-(m-1)*dk/ny):end) = (rand()*2-1)*0.1;
    end
end

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu*nu)*1.0;
Psi    = eye(N *ny)*1.0;

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(nu,kk);
y = zeros(ny,kk);
ys = zeros(ny,kk);

%% Macierze wyznaczane offline
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

% Macierz K
K = (M'*Psi*M+Lambda)^(-1)*M'*Psi;

%% Macierze dla wersji minimalistycznej algorytmu
Kyzad = zeros(nu,ny);   % r,m
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

%% Generacja macierzy
gpc2x2_matlab_to_C

%% Symulacja
du_diff = zeros(nu,kk);
for k = kp:kk
    % symulacja obiektu regulacji
    for m=1:ny
        for n=1:nu
            for i=1:nb
                if(k-i>=1)
                    ys(m,k) = ys(m,k) + b(m,n,i)*u(n,k-i);
                end
            end
        end
        for i=1:na
            if(k-i>=1)
                ys(m,k) = ys(m,k) - a(m,i)*ys(m,k-i);
            end
        end            
    end 
    
    % wprowadzanie zak³óceñ    
    y(:,k) = ys(:,k) + normrnd(mu,sigma)';
     
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
    
    u(:,k) = u(:,k-1)+du_po;
    
    for n=1:nu
        if(u(n,k)>umax); u(n,k) = umax; end
        if(u(n,k)<umin); u(n,k) = umin; end
    end
end

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
figure;
for m=1:ny
    for n=1:nu
        subplot(ny,nu,((m-1)*nu+n));
        plot(Tp*((1:size(S,3))-1),squeeze(S(m,n,:)));
        xlim([0,4]);
        if(m==1)
            ylim([0,5]);
        else
            ylim([0,2]);
        end
    end
end

figure;
plot(y'); hold on;
stairs(yzad','k--'); hold off;
title('Wartoœci wyjœciowe i zadane w czasie');

figure;
stairs(u');
title('Wartoœci sterowania w czasie');

% figure;
% stairs(du_diff');
% title('Wartoœci b³êdu w czasie');

%% Funkcje do wyznaczania minimalnej postaci algorytmu GPC
function out = fun_e(p,j,m,n)
    % wartoœci N, a, b, na, nb, nu, ny musz¹ ju¿ byæ w workspace'ie
    global N a b na nb nu ny E
    
    o = 0;
    if(isempty(E))
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
    % wartoœci N, a, b, na, nb, nu, ny musz¹ ju¿ byæ w workspace'ie
    global N a b na nb nu ny F
    
    o = 1;
    if(isempty(F))
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
    % wartoœci N, a, b, na, nb, nu, ny musz¹ ju¿ byæ w workspace'ie
    global N a b na nb nu ny G
    
    o = N;
    if(isempty(G))
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