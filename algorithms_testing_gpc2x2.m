clearvars
%close all

%% Algorytm DMC 2x2 (benchmark)
% Obiekt regulacji
%      input 1        input 2
Gs = [tf( 1,[.7 1]), tf( 5,[.3 1]);  % output 1
      tf( 1,[.5 1]), tf( 2,[.4 1])]; % output 2
Tp = 0.005;
Gz = c2d(Gs,Tp,'zoh');
ny = 2;
nu = 2;

% Based on STP script, pages 79 (tab. 3.1) and 86 (eq. 3.49)
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

% Horyzonty predykcji i sterowania
N  = 10; 
Nu = 5;

% Wartoœci trajektorii zadanej
yzad(1:ny,1:2000) =  .0;
yzad(1,   1: end) = -.8; yzad(2, 200: end) = -.1;
yzad(1, 400: end) = -.1; yzad(2, 600: end) = -.2;
yzad(1, 800: end) = -.3; yzad(2,1000: end) =  .1;
yzad(1,1200: end) =  .0; yzad(2,1400: end) =  .0;
yzad(1,1600: end) = -.2; yzad(2,1800: end) =  .2;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1+1;
kk = size(yzad,2);

% Macierze Lambda oraz Psi -- wagi funkcji kosztów (sta³e na horyzoncie
% predykcji/sterowania i równe dla ka¿dego wyjœcia/wejœcia)
Lambda = eye(Nu*nu);
Psi    = eye(N *ny);

% Wektory wartoœci sterowania oraz wyjœcia obiektu regulacji
u = zeros(nu,kk);
y = zeros(ny,kk);

% W³asnorêcznie wyznaczana odpowiedŸ skokowa
S = zeros(ny,nu,N);
for k = 1:size(S,3)
    % symulacja obiektu regulacji
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

M = cell(N,Nu);

% Matrix M
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

K = (M'*Psi*M+Lambda)^(-1)*M';
Knu = K(1:nu,:);
    
%% Symulacja
% wb = waitbar(0,'Simulation progress...');
for k = 10:kk
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
        %y(m,k) = y(m,k) + (rand()-.5)/500;
    end 
    
    % Wyznaczanie ymod
    ymod = zeros(ny,1);
    for m=1:ny
        for n=1:nu
            for i=1:nb
                if(k-i>=1)
                    ymod(m) = ymod(m) + b(m,n,i)*u(n,k-i);
                end
            end
        end
        for i=1:na
            if(k-i>=1)
                ymod(m) = ymod(m) - a(m,i)*y(m,k-i);
            end
        end            
    end 
    
    % Wyznaczanie d
    d = y(:,k)-ymod;
    
    % Wyznaczanie Y0
    Y0=zeros(ny,N);
    for m=1:ny
        for p=1:N
            Y0(m,p) = d(m);
            for n=1:nu
                for i=1:nb
                    if(-i+p<=-1)
                        Y0(m,p) = Y0(m,p) + b(m,n,i)*u(n,k-i+p);
                        fprintf('Y0(%d,%d) += b(%d,%d,%d)*u(%d,%d) = %f*%f;\n',m,p,m,n,i,n,k-i+p,b(m,n,i),u(n,k-i+p));
                    else
                        Y0(m,p) = Y0(m,p) + b(m,n,i)*u(n,k-1);
                        fprintf('Y0(%d,%d) += b(%d,%d,%d)*u(%d,%d) = %f*%f;\n',m,p,m,n,i,n,k-1,b(m,n,i),u(n,k-1));
                    end
                end
            end
            for i=1:na
                if(-i+p<=0)
                    Y0(m,p) = Y0(m,p) - a(m,i)*y(m,k-i+p);
                    fprintf('Y0(%d,%d) -= a(%d,%d)*y(%d,%d) = %f*%f;\n',m,p,m,i,m,k-i+p,a(m,i),y(m,k-i+p));
                else
                    Y0(m,p) = Y0(m,p) - a(m,i)*Y0(m,-i+p);
                    fprintf('Y0(%d,%d) -= a(%d,%d)*Y0(%d,%d) = %f*%f;\n',m,p,m,i,m,-i+p,a(m,i),Y0(m,-i+p));
                end
            end            
        end
    end
    
    Y0 = reshape(Y0,[],1);
    
    Yzad = repmat(eye(ny),N,1)*yzad(:,k); % Yzad sta³e na horyzoncie predykcji
%     Yzad = reshape(yzad(:,min(k+(1:N),kk)),[],1); % Yzad zmienne na horyzoncie predykcji
  
    du = Knu*(Yzad-Y0);
   
    u(:,k) = u(:,k-1)+du;
    for n=1:nu
        if(u(n,k)>umax); u(n,k) = umax; end
        if(u(n,k)<umin); u(n,k) = umin; end
    end
%     waitbar(k/kk,wb);
end
% close(wb)

%% Rysownie przebiegów trajektorii wyjœcia, zadanej oraz sterowania
figure;
plot(y'); hold on;
stairs(yzad','k--'); hold off;

figure;
stairs(u');