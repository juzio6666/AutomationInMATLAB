clearvars;
global Hp Hs
D=50; Hp=10; Hs=2; %wsplambda=[0.075 0.075]; wspmi=[1 1];
wsplambda=[1 1]; wspmi=[1 1];
ilwe=2; ilwy=2;

clip = [];

na=2; nb=2;
Tp=0.005;%okres pr�bkowania
K=1; T=0.7; alfa=exp(-Tp/T); B11=K*(1-alfa); A11=-alfa; % G11=1/0,7s+1
K=5; T=0.3; alfa=exp(-Tp/T); B12=K*(1-alfa); A12=-alfa; % G12=5/0,3s+1
K=1; T=0.5; alfa=exp(-Tp/T); B21=K*(1-alfa); A21=-alfa; % G21=1/0,5s+1
K=2; T=0.4; alfa=exp(-Tp/T); B22=K*(1-alfa); A22=-alfa; % G22=2/0,4s+1

as(1,1)=A11+A12; as(1,2)=A11*A12;
as(2,1)=A21+A22; as(2,2)=A21*A22;
% %b(m,n,i), m - wyj�cie, n - wej�cie, i - przesuni�cie
bs(1,1,1)=B11; bs(1,1,2)=B11*A12;
bs(1,2,1)=B12; bs(1,2,2)=B12*A11;
bs(2,1,1)=B21;
bs(2,1,2)=B21*A22;
bs(2,2,1)=B22; bs(2,2,2)=B22*A21;

% for m=1:ilwy;
%     for n=1:ilwe;
%         for i=1:nb;
%             fprintf('bs[%d][%d][%d] = %.9f;\n',m-1,n-1,i-1,bs(m,n,i));
%         end;
%     end;
%     for i=1:na;
%         fprintf('as[%d][%d] = %.9f;\n',m-1,i-1,as(m,i));
%     end;
% end;

%model=obiekt
%macierz ci�g�ych transmitancji -> macierz dyskretnych transmitancji
Tp=0.0001;%okres pr�bkowania
K=1; T=0.7; alfa=exp(-Tp/T); B11=K*(1-alfa); A11=-alfa; % G11=1/0,7s+1
K=5; T=0.3; alfa=exp(-Tp/T); B12=K*(1-alfa); A12=-alfa; % G12=5/0,3s+1
K=1; T=0.5; alfa=exp(-Tp/T); B21=K*(1-alfa); A21=-alfa; % G21=1/0,5s+1
K=2; T=0.4; alfa=exp(-Tp/T); B22=K*(1-alfa); A22=-alfa; % G22=2/0,4s+1

%sprowadzenie do wsp�lnego mianownika (r�wnania r�nicowe)
na=2; nb=2;
% %a(m,i), m - wyj�cie, i - przesuni�cie
a(1,1)=A11+A12; a(1,2)=A11*A12;
a(2,1)=A21+A22; a(2,2)=A21*A22;
% %b(m,n,i), m - wyj�cie, n - wej�cie, i - przesuni�cie
b(1,1,1)=B11; b(1,1,2)=B11*A12;
b(1,2,1)=B12; b(1,2,2)=B12*A11;
b(2,1,1)=B21;
b(2,1,2)=B21*A22;
b(2,2,1)=B22; b(2,2,2)=B22*A21;

clip = [clip sprintf('#define HP %d\n',Hp)];
clip = [clip sprintf('#define HS %d\n',Hs)];

clip = [clip sprintf('#define ilwy %d\n',ilwy)];
clip = [clip sprintf('#define ilwe %d\n',ilwe)];

clip = [clip sprintf('#define na %d\n',na)];
clip = [clip sprintf('#define nb %d\n',nb)];

clip = [clip sprintf('float b[ilwy][ilwe][nb];\n')];
clip = [clip sprintf('float a[ilwy][na];\n')];
clip = [clip sprintf('void Model_Init(void){\n')];
for m=1:ilwy;
    for n=1:ilwe;
        for i=1:nb;
            clip = [clip sprintf('\tb[%d][%d][%d] = %.9f;\n',m-1,n-1,i-1,b(m,n,i))];
        end;
    end;
    for i=1:na;
        clip = [clip sprintf('\ta[%d][%d] = %.9f;\n',m-1,i-1,a(m,i))];
    end;
end;
clip = [clip sprintf('}\n')];

kp=max(na,nb)+1;%pocz�tek symulacji
%warunki pocz�tkowe
u(1:2,1:kp-1)=0; y(1:2,1:kp-1)=0;

%trajektoria jak w ksi��ce
yzad(1:2,1:9)=0;
yzad(1,kp:200)=0;
yzad(2,kp:200)=0;
yzad(1,201:1000)=100;
yzad(2,201:1000)=  0;
yzad(1,1001:2000)=-100;
yzad(2,1001:2000)=  0;

[q, kk]=size(yzad);

%odpowied� skokowa
%s(1,ilwy,1,ilwe,1,Hp);
s = zeros(ilwy,ilwe,Hp);
for k=1:D
    i1=min(k,nb); i2=min(k-1,na);
    for m=1:ilwy;
        for n=1:ilwe
            s(m,n,k)=0.0;
            for i=1:i1;
                s(m,n,k)=s(m,n,k)+b(m,n,i);
            end;
            for i=1:i2;
                s(m,n,k)=s(m,n,k)-a(m,i)*s(m,n,k-i);
            end;
        end;
    end;
end;

for m=1:ilwy;
    for n=1:ilwe
        s(m,n,D+1:D+Hp)=s(m,n,D);
    end;
end;
%
G=zeros(ilwy*Hp,ilwe*Hs);
%pierwsza kolumna z�o�ona z macierzy s
i=1;%wiersz
j=1;%kolumna
for k=1:Hp;
    for m=1:ilwy;
        for n=1:ilwe;
            j=n;
            G(i,j)=s(m,n,k);
        end;
        i=i+1;
    end;
end;
%kolejne kolumny
for p=1:Hs-1;
    for i=p*ilwy+1:ilwy*Hp;
        for j=p*ilwe+1:p*ilwe+ilwe;
            G(i,j)=G(i-p*ilwy,j-p*ilwe);
        end;
    end;
end;

MP=zeros(ilwy*Hp,ilwe*(D-1));
%wiersz=1; kolumna=1;
for i=1:Hp;
    for j=1:D-1;
        for m=1:ilwy;
            for n=1:ilwe
                wiersz=2*i-1+m-1; kolumna=2*j-1+n-1;
                MP(wiersz,kolumna)=s(m,n,i+j)-s(m,n,j);
            end;
        end;
    end;
end;

%macierze MI LAMBDA
MI = zeros(Hp*ilwy);
i=1;
for p=1:Hp;
    for m=1:ilwy;
        MI(i,i)=wspmi(m); i=i+1;
    end;
end;
i=1;
LAMBDA = zeros(Hs*ilwe);
for p=1:Hs;
    for n=1:ilwe;
        LAMBDA(i,i)=wsplambda(n); i=i+1;
    end;
end;

K=(G'*MI*G+LAMBDA)\G'*MI;
tmpK = K';
tmpTxt = 'float K[(HS*ilwe)*(HP*ilwy)] = {';
clip = [clip sprintf(tmpTxt)];
for x1=1:size(K,1)
    if x1~=1
        clip = [clip sprintf('%s',ones(size(tmpTxt))*' ')];
    end
    for x2=1:size(K,2)
        clip = [clip sprintf('%.9f',K(x1,x2))];
        if x2==size(K,2) && x1==size(K,1)
            clip = [clip sprintf('};')];
        else 
            clip = [clip sprintf(',')];
        end
    end
    clip = [clip sprintf('\n')];
end

DUP=zeros(ilwe*(D-1),1);
Yzad=zeros(ilwy*Hp,1);
Yk=zeros(ilwy*Hp,1);
Y0=zeros(ilwy*Hp,1);
y0=zeros(ilwy,kk);
for k=kp:kk;
    %symulacja obiektu
    for m=1:ilwy;
        y(m,k)=0.0;
        for n=1:ilwe;
            for i=1:nb;
                y(m,k)=y(m,k)+b(m,n,i)*u(n,k-i);
            end;
        end;
        for i=1:na;
            y(m,k)=y(m,k)-a(m,i)*y(m,k-i);
        end;
    end;
    
    % %tu mo�na doda� ewentualne niemierzalne zak��cenia
%     if k>=30
%         y(1,k)=y(1,k)+0.01;
%         y(2,k)=y(2,k)-0.01;
%     end;
    
    %algorytm gpc
    for m=1:ilwy;
        ymod=0.0;
        for n=1:ilwe;
            for i=1:nb;
                ymod=ymod+b(m,n,i)*u(n,k-i);
            end;
        end;
        for i=1:na;
            ymod=ymod-a(m,i)*y(m,k-i);
        end;
        d=y(m,k)-ymod;
        for p=1:Hp;
            y0(m,p)=d;
            for n=1:ilwe;
                for i=1:nb;
                    if (-i+p>=0)
                        liczba=u(n,k-1);
                    else
                        liczba=u(n,k-i+p);
                    end;
                    y0(m,p)=y0(m,p)+b(m,n,i)*liczba;
                end;
            end;
            for i=1:na;
                if (-i+p>=1)
                    liczba=y0(m,p-i);
                else
                    liczba=y(m,k+p-i);
                end;
                y0(m,p)=y0(m,p)-a(m,i)*liczba;
            end;
        end;%p=1...Hp
    end;%m=1...ilwy

    l=1;
    for p=1:Hp;
        Y0(l)=y0(1,p); l=l+1;
        Y0(l)=y0(2,p); l=l+1;
    end;
    
    l=1;
    for p=1:Hp;
        for m=1:ilwy
            Yzad(l)=yzad(m,k);
            l=l+1;
        end;
    end;
    
    DU=K*(Yzad-Y0);
    
    u(1,k)=DU(1)+u(1,k-1);
    u(2,k)=DU(2)+u(2,k-1);
end;
% txt=sprintf('GPC: Hp=%ld, Hs=%ld',Hp,Hs);

% figure;
% subplot(2,1,1); stairs(u(1,1:kk)); grid; xlabel('k'); ylabel('u_1'); title(txt);
% subplot(2,1,2); stairs(u(2,1:kk)); grid; xlabel('k'); ylabel('u_2');
% figure;
% subplot(2,1,1); stairs(yzad(1,1:kk),'r'); hold on; plot(y(1,1:kk),'b'); grid; xlabel('k'); ylabel('y_1^{zad}, y_1'); title(txt);
% subplot(2,1,2); stairs(yzad(2,1:kk),'r'); hold on; plot(y(2,1:kk),'b'); grid; xlabel('k'); ylabel('y_2^{zad}, y_2');
fprintf(clip);
clipboard('copy',clip);