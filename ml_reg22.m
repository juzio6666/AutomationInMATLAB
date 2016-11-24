clear all;
D=50; Hp=10; Hs=5; %wsplambda=[0.075 0.075]; wspmi=[1 1];
wsplambda=[0.075 0.15]; wspmi=[1 2];

typalgorytmu=2;%1 - DMC, 2 - GPC
optymalizacja=1;%0 - wersja analityczna, 1 - wersja numeryczna
u1min=-2.5-0*5; u2min=-0.6-0*5;
u1max=2.5+0*5; u2max=0.6+0*5;
du1max=0.6; du2max=0.1; jestogrdumax=1;

ilwe=2; ilwy=2;

%model=obiekt

%macierz ci¹g³ych transmitancji -> macierz dyskretnych transmitancji
% Tp=0.03;%okres próbkowania
% %G11=1/0,7s+1
% K=1; T=0.7;
% alfa=exp(-Tp/T); B11=K*(1-alfa); A11=-alfa;
% %G12=5/0,3s+1
% K=5; T=0.3;
% alfa=exp(-Tp/T); B12=K*(1-alfa); A12=-alfa;
% %G21=1/0,5s+1
% K=1; T=0.5;
% alfa=exp(-Tp/T); B21=K*(1-alfa); A21=-alfa;
% %G22=2/0,4s+1
% K=2; T=0.4;
% alfa=exp(-Tp/T); B22=K*(1-alfa); A22=-alfa;

%sprowadzenie do wspólnego mianownika (równania ró¿nicowe)
na=2; nb=2;
% %a(m,i), m - wyjœcie, i - przesuniêcie
% a(1,1)=A11+A12; a(1,2)=A11*A12;
% a(2,1)=A21+A22; a(2,2)=A21*A22;
% %b(m,n,i), m - wyjœcie, n - wejœcie, i - przesuniêcie
% b(1,1,1)=B11; b(1,1,2)=B11*A12;
% b(1,2,1)=B12; b(1,2,2)=B12*A11;
% b(2,1,1)=B21;
% b(2,1,2)=B21*A22;
% b(2,2,1)=B22; b(2,2,2)=B22*A21;
%wartoœci liczbowe
a(1,1)=-1.862885; a(1,2)=0.866877;%y1
a(2,1)=-1.869508; a(2,2)=0.873715;%y2
b(1,1,1)=0.041951; b(1,1,2)=-0.037959; b(1,2,1)=0.475812; b(1,2,2)=-0.455851;
b(2,1,1)=0.058235; b(2,1,2)=-0.054027; b(2,2,1)=0.144513; b(2,2,2)=-0.136097;

kp=max(na,nb)+1;%pocz¹tek symulacji

%warunki pocz¹tkowe
u(1:2,1:kp-1)=0; y(1:2,1:kp-1)=0;

%pojedynczy skok
% yzad(1:2,1:2)=0;
% yzad(1,3:150)=1;
% yzad(2,3:150)=2;

%trajektoria jak w ksi¹¿ce
yzad(1:2,1:9)=0;
yzad(1,10:49)=0.7;
yzad(1,50:150)=-0.7;
yzad(2,10:99)=-0.5;
yzad(2,100:150)=0.5;

[q kk]=size(yzad);

%odpowiedŸ skokowa
%s(1,ilwy,1,ilwe,1,Hp);
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
%G=zeros(ilwy*Hp,ilwe*Hs);
%pierwsza kolumna z³o¿ona z macierzy s
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
                %disp(sprintf('i=%d, j=%d, m=%d, n=%d, MP(%d,%d)',i,j,m,n,wiersz,kolumna));
                MP(wiersz,kolumna)=s(m,n,i+j)-s(m,n,j);
            end;
        end;
    end;
end;
%MP(i,j)=s(i+j)-s(j);%1x1

%macierze MI LAMBDA
i=1;
for p=1:Hp;
    for m=1:ilwy;
        MI(i,i)=wspmi(m); i=i+1;
    end;
end;
i=1;
for p=1:Hs;
    for n=1:ilwe;
        LAMBDA(i,i)=wsplambda(n); i=i+1;
    end;
end;

K=inv(G'*MI*G+LAMBDA)*G'*MI;

DUP=zeros(ilwe*(D-1),1);
Yzad=zeros(ilwy*Hp,1);
Yk=zeros(ilwy*Hp,1);
Y0=zeros(ilwy*Hp,1);

H_qp=2*(G'*MI*G+LAMBDA);
H_qp=0.5*(H_qp+H_qp');
opcje_qp=optimset('Algorithm','active-set','LargeScale','off','Display','off');
J=zeros(2*Hs,2*Hs);
w=1; l=1;
for i=1:Hs;
    k=1;
    for j=1:l;
        J(w:w+ilwe-1,k:k+ilwe-1)=eye(ilwe,ilwe);
        k=k+ilwe;
    end;
    w=w+ilwe; l=l+1;
end;
Aogr_qp=[-J; J];
bogr_qp=zeros(2*Hs,1);

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
    
    % %tu mo¿na dodaæ ewentualne niemierzalne zak³ócenia
%     if k>=30
%         y(1,k)=y(1,k)+0.01;
%         y(2,k)=y(2,k)-0.01;
%     end;
    
    %algorytm dmc
    if typalgorytmu==1
        l=1;
        for p=1:D-1;
            for n=1:ilwe
                DUP(l)=0;
                if k-p>0
                    DUP(l)=DUP(l)+u(n,k-p);
                end;
                if k-p-1>0
                    DUP(l)=DUP(l)-u(n,k-p-1);
                end;
                l=l+1;
            end;
            
        end;
        
        l=1;
        for p=1:Hp;
            for m=1:ilwy
                Yk(l)=y(m,k);
                l=l+1;
            end;
        end;
        
        Y0=Yk+MP*DUP;
    end;
    
    %algorytm gpc
    if typalgorytmu==2
        for m=1:ilwy;
            ymod=0.0;
            for n=1:ilwe;
                for i=1:na;
                    ymod=ymod+b(m,n,i)*u(n,k-i);
                end;
            end;
            for i=1:nb;
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
    end;
    
    l=1;
    for p=1:Hp;
        for m=1:ilwy
            Yzad(l)=yzad(m,k);
            l=l+1;
        end;
    end;
    
    if optymalizacja==0
        DU=K*(Yzad-Y0);
    end;
    if optymalizacja==1
        f_qp=-2*G'*MI*(Yzad-Y0);
        l=1;
        for p=0:Hs-1;
            bogr_qp(l)=-u1min+u(1,k-1); l=l+1;
            bogr_qp(l)=-u2min+u(2,k-1); l=l+1;
        end;
        for p=0:Hs-1;
            bogr_qp(l)=u1max-u(1,k-1); l=l+1;
            bogr_qp(l)=u2max-u(2,k-1); l=l+1;
        end;
        if jestogrdumax==0
            vlb_qp=[];
            vub_qp=[];
        end;
        if jestogrdumax==1
            l=1;
            for p=0:Hs-1;
                vlb_qp(l)=-du1max; vub_qp(l)=du1max; l=l+1;
                vlb_qp(l)=-du2max; vub_qp(l)=du2max; l=l+1;
            end;
        end;
        %du=qp(H_qp,f_qp,Aogr_qp,bogr_qp,vlb_qp,vub_qp);
        %X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS)
        DU=quadprog(H_qp,f_qp,Aogr_qp,bogr_qp,[],[],vlb_qp,vub_qp,[],opcje_qp);
        
        if k>=10 & 0
            k
            H_qp
            f_qp
            Aogr_qp
            bogr_qp
            %Aogr_qp*zeros(2*Hs,1)-bogr_qp
            pause
        end;
        
    end;
    
    u(1,k)=DU(1)+u(1,k-1);
    u(2,k)=DU(2)+u(2,k-1);
    
    %txt=sprintf('k=%d',k); disp(txt);
    %txt=sprintf(' y10(k)=%e, %e, %e',y0(1,1),y0(1,2),y0(1,3)); disp(txt);    
    %txt=sprintf(' y20(k)=%e, %e, %e',y0(2,1),y0(2,2),y0(2,3)); disp(txt);        
    %txt=sprintf(' u1(k)=%e, u1(k)=%e',u(1,k),u(2,k)); disp(txt);    
    %pause
end;
e1=(yzad(1,1:kk)-y(1,1:kk))*(yzad(1,1:kk)-y(1,1:kk))';
e2=(yzad(2,1:kk)-y(2,1:kk))*(yzad(2,1:kk)-y(2,1:kk))';
if typalgorytmu==1
    txt=sprintf('DMC: optymalizacja=%d, D=%ld, Hp=%ld, Hs=%ld, e1=%e, e2=%e',optymalizacja,D,Hp,Hs,e1,e2);
end;
if typalgorytmu==2
    txt=sprintf('GPC: optymalizacja=%d, Hp=%ld, Hs=%ld, e1=%e, e2=%e',optymalizacja,Hp,Hs,e1,e2);
end;

figure;
subplot(2,1,1); stairs(u(1,1:kk)); grid; xlabel('k'); ylabel('u_1'); title(txt);
if optymalizacja==1
    axis([0 150 u1min u1max]);
end;
subplot(2,1,2); stairs(u(2,1:kk)); grid; xlabel('k'); ylabel('u_2');
if optymalizacja==1
    axis([0 150 u2min u2max]);
end;
figure;
subplot(2,1,1); stairs(yzad(1,1:kk),'r'); hold on; plot(y(1,1:kk),'b'); grid; xlabel('k'); ylabel('y_1^{zad}, y_1'); title(txt);
subplot(2,1,2); stairs(yzad(2,1:kk),'r'); hold on; plot(y(2,1:kk),'b'); grid; xlabel('k'); ylabel('y_2^{zad}, y_2');



