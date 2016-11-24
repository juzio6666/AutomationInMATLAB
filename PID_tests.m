function [u,y,z,awup,test]=PID_tests(T,K,Ti,Td,Tv)
[a,b,umin,umax,~,~]= obiekt3();
na=length(a); nb=length(b); 
kp=max(na,nb)+1; kp=5; kk=200;
u(1:kp-1)=0; u(kp:kk)=0; 
uw(1:kp-1)=0; uw(kp:kk)=0; 
y(1:kp-1)=0; y(kp:kk)=0; 
e(1:kp-1)=0; e(kp:kk)=0; 
z(1:kp-1)=0; z(kp:kk)=0; % yzad
up(1:kk)=0;
ud(1:kk)=0;
ui(1:kk)=0;
u1(1:kk)=0;
u2(1:kk)=0;
awup(1:kk)=0;
z(10:kk) = 1;
test(1:kk)=0;

for k=kp:kk;
    %symulacja obiektu
    y(k)=0;
    %fprintf('y(k) = ');
    for i=1:nb
        y(k)=y(k)+b(i)*uw(k-i);
    %    fprintf('+%f*u(k-%d)',b(i),i);
    end;
    for i=1:na
        y(k)=y(k)-a(i)*y(k-i);
    %    fprintf('-%f*y(k-%d)',a(i),i);
    end;
    %fprintf('\n');
    
    e(k)=z(k)-y(k);
    %e(k)=0.1;          % do wyznaczania czasu zdwojenia (PI)
    %e(k)=e(k-1)+0.1;   % do wyznaczania czasu wyprzedzenia (PD)
    
    % Podejscie pierwsze do wyznaczenia wartosci sterowania
    Tibis = 1/(1/Ti + awup(k-1)/Tv);
    test(k)=Tibis;
    r2 = K*Td/T;                    
    r1 = K*(T/(2*Tibis) - 2*Td/T - 1); 
    r0 = K*(1 + T/(2*Tibis) + Td/T);   
    u1(k) = u(k-1)+r2*e(k-2)+r1*e(k-1)+r0*e(k);
    
    % Podejscie drugie do wyznaczenia wartosci sterowania
    up(k) = K*e(k);
    ud(k) = K*Td*(e(k)-e(k-1))/T;
    ui(k) = ui(k-1)+K*(1/Ti+awup(k-1)/Tv)*T*(e(k-1)+e(k))/2;
    u2(k) = up(k)+ud(k)+ui(k);
    
    u(k) = u2(k);
    uw(k) = max(min(u(k),umax),umin);
    awup(k) = uw(k)-u(k);
    u(k) = uw(k);
end;
end

function [a,b,umin,umax,dumax,Tp] = obiekt1()
    Tp = 0.03;
    K=3; T=30;
    [ld,md]=c2dm(K,[T 1],Tp,'zoh');
    b(1)=ld(2);
    a(1)=md(2);
    umin=-0.5; 
    umax=0.5;
    dumax=0.1;
end

function [a,b,umin,umax,dumax,Tp] = obiekt2()
    Tp = 0.03;
    K=2; 
    T1=0.03; 
    T2=0.5;
    [ld,md]=c2dm(K,[T1*T2 T1+T2 1],Tp,'zoh');
    b(1)=ld(2);
    b(2)=ld(3);
    a(1)=md(2);
    a(2)=md(3);
    %umin=-1.0;
    %umax=1.0;
    umin=-Inf;
    umax=Inf;
    dumax=0.01;
end

function [a,b,umin,umax,dumax,Tp] = obiekt3()
    Tp = 2;
    a = [-1.6375, 0.67003];
    b = [0, 0, 0.035, 0.0307];
    umin = -Inf;
    umax = Inf;
    dumax = Inf;
end