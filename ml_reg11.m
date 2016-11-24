clear all;
D=50; Hp=10; Hs=5; lambda=2.5;

optymalizacja=1;%0 - wersja analityczna, 1 - wersja numeryczna
typobiektu=2;%1 lub 2 inercje
Tp=1;%okres próbkowania

if typobiektu==1;
    a(1)=-0.904837; b(1)=0.285487;
    
    umin=-0.5; umax=0.5;
    dumax=0.1; jestogrdumax=1;        
end;

if typobiektu==2;
    a(1)=-1.621368; a(2)=0.648344;
    b(1)=0.028919; b(2)=0.025031;
    
    umin=-0.7; umax=0.7;
    dumax=0.2; jestogrdumax=1;    
end;
na=length(a); nb=length(b); kp=max(na,nb)+1;

yzad(1:kp-1)=0; yzad(kp:100)=1;
kk=length(yzad);

%odpowiedŸ skokowa
for k=1:D;
    s(k)=0;
    for i=1:min(k,nb);
        s(k)=s(k)+b(i);
    end;
    for i=1:min(k-1,na);
        s(k)=s(k)-a(i)*s(k-i);
    end;
end;
s(D+1:D+Hp)=s(D);

G=zeros(Hp,Hs);
for i=1:Hp;
    G(i,1)=s(i);
end;
for i=2:Hs;
    G(i:Hp,i)=G(1:Hp-i+1,1);
end;

MP=zeros(Hp,D-1);
for i=1:Hp;
    for j=1:D-1;
        MP(i,j)=s(i+j)-s(j);
    end;
end;

K=inv(G'*G+lambda*eye(Hs,Hs))*G';      
DUP=zeros(D-1,1);
Yzad=zeros(Hp,1);
Yk=zeros(Hp,1);
Y0=zeros(Hp,1);

H_qp=2*(G'*G+lambda*eye(Hs,Hs));
H_qp=0.5*(H_qp+H_qp');
opcje_qp=optimset('Algorithm','interior-point-convex','LargeScale','off','Display','off'); 
J=zeros(Hs,Hs);
for i=1:Hs;
    J(i,1:i)=1;
end;
Aogr_qp=[-J; J];
bogr_qp=zeros(2*Hs,1);

%warunki pocz¹tkowe
u(1:kp-1)=0; y(1:kp-1)=0;
for k=kp:kk;
    %symulacja obiektu
    y(k)=0;
    for i=1:nb
        y(k)=y(k)+b(i)*u(k-i);
    end;
    for i=1:na
        y(k)=y(k)-a(i)*y(k-i);
    end;
    
    %niemierzalne zak³ócenie
%      if k>=40
%          y(k)=y(k)+0.01;
%      end;
    
    %algorytm gpc
    ym=0;
    for i=1:nb;
        ym=ym+b(i)*u(k-i);
    end;
    for i=1:na;
        ym=ym-a(i)*y(k-i);
    end;
    ddmc=y(k)-ym;

    %odpowiedŸ swobodna
    for p=1:Hp;
        y0(p)=ddmc;
        for i=1:nb;
            if -i+p>=0
                y0(p)=y0(p)+b(i)*u(k-1);
            end;
            if -i+p<0
                y0(p)=y0(p)+b(i)*u(k-i+p);
            end;
        end;
        for i=1:na;
            if -i+p>=1
                y0(p)=y0(p)-a(i)*y0(-i+p);
            end;
            if -i+p<1
                y0(p)=y0(p)-a(i)*y(k-i+p);
            end;
        end;
    end;
    Y0=y0';
        
    Yzad=yzad(k)*ones(Hp,1);
    if optymalizacja==0
        DU=K*(Yzad-Y0);
    end;
    if optymalizacja==1
        f_qp=-2*G'*(Yzad-Y0);
        l=1;
        for p=0:Hs-1;
            bogr_qp(l)=-umin+u(k-1); l=l+1;
        end;
        for p=0:Hs-1;
            bogr_qp(l)=umax-u(k-1); l=l+1;
        end;
        if jestogrdumax==0
            vlb_qp=[];
            vub_qp=[];
        end;
        if jestogrdumax==1
            vlb_qp=-dumax*ones(Hs,1);
            vub_qp=dumax*ones(Hs,1);
        end;
        %du=qp(H_qp,f_qp,Aogr_qp,bogr_qp,vlb_qp,vub_qp);
        %X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS)
        DU=quadprog(H_qp,f_qp,Aogr_qp,bogr_qp,[],[],vlb_qp,vub_qp,[],opcje_qp);
        %[k DU']
        %f_qp
        %pause
        %rozwi¹zanie pocz¹tkowe zerowe - spe³nia ograniczenia
%         Aogr_qp*zeros(Hs,1)-bogr_qp
%         bogr_qp
%         pause
    end;
    u(k)=DU(1)+u(k-1);
end;

%wyniki symulacji
e=(yzad(1:kk)-y(1:kk))*(yzad(1:kk)-y(1:kk))';
txt=sprintf('GPC: optymalizacja=%d, Hp=%ld, Hs=%ld, e=%e',optymalizacja,Hp,Hs,e);

figure; stairs(u,'b');
xlabel('k');
ylabel('u');
title(txt);
grid;

figure; stairs(yzad,'r');
hold on; plot(y,'b');
xlabel('k');
ylabel('y_{zad}, y');
title(txt);
grid;