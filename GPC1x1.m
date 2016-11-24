function [y,u,yzad]=GPC1x1()
    D=50; Hp=10; Hs=5; lambda=2.5;
    typobiektu=2;%1 lub 2 inercje

    if typobiektu==1;
        a(1)=-0.904837; b(1)=0.285487;
    
        umin=-0.5; umax=0.5;
        dumax=0.1; jestogrdumax=0;        
    end;

    if typobiektu==2;
        a(1)=-1.621368; a(2)=0.648344;
        b(1)=0.028919; b(2)=0.025031;
    
        umin=-0.7; umax=0.7;
        dumax=0.2; jestogrdumax=0;
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

    K=inv(G'*G+lambda*eye(Hs,Hs))*G';      
    DUP=zeros(D-1,1);
    Yzad=zeros(Hp,1);
    Yk=zeros(Hp,1);
    Y0=zeros(Hp,1);

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
        DU=K*(Yzad-Y0);
        u(k)=DU(1)+u(k-1);
    end;
end