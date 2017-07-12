clearvars -except fans
global N a b na nb u y kp kk
global Y0temp Y0temp2
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

t1=1.0; t2=2.0;
alfa1=exp(-1.0/t1); alfa2=exp(-1.0/t2);
a(1)=-alfa1-alfa2; a(2)=alfa1*alfa2;
b(1)=(t1/(t1-t2))*(1.0-alfa1)+(t2/(t2-t1))*(1.0-alfa2);
b(2)=(t1/(t1-t2))*(alfa1-1.0)*alfa2+(t2/(t2-t1))*(alfa2-1.0)*alfa1;
na=2; nb=2;

% Ograniczenia
umax =  1;
umin = -1;

% Horyzonty predykcji i sterowania
N = 5; 
Nu = 2;

% Wartoœci trajektorii zadanej
yzad(   1:2000) =  .0;
yzad( 200: end) = -.1;
yzad( 400: end) =  .1;
yzad( 600: end) =  .2;
yzad( 800: end) =  .3;
yzad(1000: end) = -.2;
yzad(1200: end) = -.3;
yzad(1400: end) =  .0;
yzad(1600: end) =  .4;
yzad(1800: end) = -.4;

% Pocz¹tkowa i koñcowa chwila symulacji
kp = max(na,nb)+1;
kk = length(yzad);

% Macierze Lambda oraz Psi -- wagi funkcji kosztów
Lambda = eye(Nu)*0.1;
Psi = eye(N)*1.0;

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
Kyzad = sum(Knu);

Ku=zeros(1,nb);
for i=1:nb
    Ku(i) = 0;
    for p=1:N
        Ku(i) = Ku(i)+K(1,p)*fun_e(p,i);
    end
end

Ky=zeros(1,na+1);
for i=0
    for p=1:N
        Ky(1) = Ky(1)+K(1,p)*fun_f(p,0);
    end
end
for i=1:na
    for p=1:N
        Ky(i+1) = Ky(i+1)+K(1,p)*fun_f(p,i);
    end
end

%% Symulacja
for k = kp:kk
    y(k) = -a*y(k-(1:length(a)))+b*u(k-(1:length(b))); % symulacja obiektu regulacji
    
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
    
%     du = Knu*(Yzad-Y0);
%     Y0temp(k,:) = -Knu*Y0;
%     Y0temp2(k,:) = - Ku*u(k-(1:nb)) - Ky*y(k-(0:na));
%     Ku=fans(1:2); Ky=fans(3:5);
    du = Kyzad*yzad(k) - Ku*u(k-(1:nb)) - Ky*y(k-(0:na));
    
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


function out = fun_e(p,j)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb 
    if(j==1 && nb==1)
        out = 0;
    elseif(j==1 && nb>1)
        out = fun_g(p,j);
    elseif(j>=2 && j<=(nb-1) && j<nb && nb>1)
        out = fun_g(p,j) - fun_g(p,j-1);
    elseif(j==nb && nb>1)
        out = -fun_g(p,j-1);
    else
        error('Coœ posz³o nie tak! Parametry wywo³ania:\np=%d; j=%d; nb=%d;\n',p,j,nb);
    end
end

function out = fun_f(p,j)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb 
    if(p == 1)
        if(j==0)
            out = 1 - a(1);
        elseif(j>=1 && j<=(na-1))
            out = a(j)-a(j+1);
        elseif(j==na)
            out = a(j);
        else
            error('Coœ posz³o nie tak1! Parametry wywo³ania:\np=%d; j=%d; na=%d;\n',p,j,na);
        end
    elseif(p>=2 && p<=N)
        if(j>=0 && j<=(na-1))
            out = fun_f(p-1,0)*fun_f(1,j)+fun_f(p-1,j+1);
        elseif(j==na)
            out = fun_f(p-1,0)*fun_f(1,j);
        else
            error('Coœ posz³o nie tak2! Parametry wywo³ania:\np=%d; j=%d; na=%d;\n',p,j,na);
        end
    else
        error('Coœ posz³o nie tak3! Parametry wywo³ania:\np=%d; j=%d; na=%d;\n',p,j,na);
    end
end

function out = fun_g(p,j)
    % wartoœci N, a, b, na, nb musz¹ ju¿ byæ w workspace'ie
    global N a b na nb 
    if(p == 1)
        if(j>=0 && j<=(nb-1))
            out = b(j+1);
        else
            error('Coœ posz³o nie tak1! Parametry wywo³ania:\np=%d; j=%d; nb=%d;\n',p,j,nb);
        end
    elseif(p>=2 && p<=N)
        if(j>=(1-p) && j<=(-1))
            out = fun_g(p-1,j+1);
        elseif(j>=0 && j<=(nb-2))
            out = fun_f(p-1,0)*fun_g(1,j)+fun_g(p-1,j+1);
        elseif(j==(nb-1))
            out = fun_f(p-1,0)*fun_g(1,j);
        else
            error('Coœ posz³o nie tak2! Parametry wywo³ania:\np=%d; j=%d; nb=%d;\n',p,j,nb);
        end
    else
        error('Coœ posz³o nie tak3! Parametry wywo³ania:\np=%d; j=%d; nb=%d;\n',p,j,nb);
    end
end