%% Algorytm DMC 2x2 (benchmark)
clear all

close all

mu   =[-0.165315345000003,-0.065108360000001]*0;
sigma=[ 0.001148139448443, 0.001082560427385]*0;

obiekt_losowy = 0;
ny = 2;
nu = 2;
%% Obiekt regulacji
if(obiekt_losowy == 0)
    inercje = 1;

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

% OdpowiedŸ skokowa 
S = step(Gz(:,:)); % wymiary D x ny x nu
S = shiftdim(S(1:end,:,:),1); % usuwanie pierwszego elementu odpowiedzi 
                              % skokowej -- KONIECZNE dla obiektu bez 
                              % opóŸnienia (dla obiektu z opóŸnieniem 
                              % pierwsze zero jest istotne)! zmieniona 
                              % zosta³a kolejnoœæ indeksów
                              
D = size(S,3);
D = 200;
%D = min(D,200); % nadpisujê ¿eby zmniejszyæ liczbê obliczeñ

% % W³asnorêcznie wyznaczana odpowiedŸ skokowa
% S = zeros(ny,nu,1000);
% for k = 1:size(S,3)
%     % symulacja obiektu regulacji
%     for m=1:ny
%         for n=1:nu
%             for i=1:min(k,nb)
%                 S(m,n,k) = S(m,n,k) + b(m,n,i)*1;
%             end
%             for i=1:min(k-1,na)
%                 S(m,n,k) = S(m,n,k) - a(m,i)*S(m,n,k-i);
%             end         
%         end   
%     end 
% end

%% Ogólne parametry algorytmu
% Horyzonty predykcji i sterowania
N  = D; 
Nu = D;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1+1;
kk = 2000;
dk = 200;

% Wartoœci trajektorii zadanej
yzad = zeros(ny,kk);
% for k=dk:dk:kk
%     for m=1:ny
%         yzad(m,(k-(m-1)*dk/ny):end) = (rand()*2-1)*0.1;
%     end
% end
yzad(1, 100:end) = -.1;
yzad(1, 500:end) =  .1;
yzad(1, 900:end) = -.3;
yzad(1,1300:end) =  .2;
yzad(1,1700:end) =  .0;

yzad(2, 300:end) = -.1;
yzad(2, 700:end) = -.2;
yzad(2,1100:end) =  .1;
yzad(2,1500:end) = -.2;
yzad(2,1900:end) =  .0;   

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu*nu)*1.0;
Psi    = eye(N *ny)*1.0;

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(nu,kk);
y = zeros(ny,kk);
ys = zeros(ny,kk);

%% Macierze wyznaczane offline
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

% Macierz Mp
Mp = cell(N,(D-1));
for row = 1:N
   for col = 1:(D-1)
        Mp{row,col} = S(:,:,min(row+col,D)) - S(:,:,col);
   end
end
Mp = cell2mat(Mp);

% Macierz K
K = (M'*Psi*M+Lambda)^(-1)*M'*Psi;

%% Macierze dla wersji minimalistycznej algorytmu
Ku=K(1:nu,:)*Mp;

Ke=zeros(nu,1);
for n=1:nu
    for m=1:ny
        Ke(n,m) = sum(K(n,m:ny:end));
    end
end

%% Generacja macierzy
dmc2x2_matlab_to_C

%% Symulacja
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
    
    % wyznaczanie dUp
    dUpp = zeros(nu,D-1);
    for p = 1:(D-1)
        if(k-p > 0); dUpp(:,p) = u(:,k-p); end
        if(k-p-1 > 0); dUpp(:,p) = dUpp(:,p)-u(:,k-p-1); end
    end
    dUp = reshape(dUpp,[],1);
   
    % wyznaczanie Yzad (sta³e na horyzoncie predykcji)
    Yzad = repmat(eye(ny),N,1)*yzad(:,k); 
   
    % wyznaczanie Y
    Y = repmat(eye(ny),N,1)*y(:,k); 
    
    
    %% wyznaczenie du (klasycznie)
    Y0 = Y+Mp*dUp;
    dU = K*(Yzad-Y0);
    du_no = dU(1:nu);
    
    %% wyznaczenie du (optymalnie)
    du = Ke*(yzad(:,k)-y(:,k)) - Ku*dUp;
    
    du_diff(1:nu,k) = du-du_no;
    
    u(:,k) = u(:,k-1)+du;
    for n=1:nu
        if(u(n,k)>umax); u(n,k) = umax; end
        if(u(n,k)<umin); u(n,k) = umin; end
    end
end

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
% figure;
% for m=1:ny
%     subplot(ny,1,m)
%     plot(y(m,:)','b'); hold on;
%     stairs(yzad(m,:)','k--'); 
%     stairs(u');
%     hold off;
% end
% title('Wartoœci wyjœciowe i zadane w czasie');

figure;
plot(y'); hold on;
stairs(yzad','k--'); hold off;
title('Wartoœci wyjœciowe i zadane w czasie');

figure;
stairs(u');
title('Wartoœci sterowania w czasie');

csvwrite('../LaTeX/DUNNO_Pomiar_czasu_algorytmow_regulacji/dane/dmc2x2sim.csv',[((1:2000)'-1)*Tp y' yzad' u']);

% figure;
% stairs(du_diff');
% title('Wartoœci b³êdu w czasie');