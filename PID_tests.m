function PID_tests()

[a,b,umin,umax,~,~]= obiekt3();
na=length(a); nb=length(b); 
kp=max(na,nb)+1; kp=5; kk=200;
u(1:kp-1)=0; u(kp:kk)=0; 
y(1:kp-1)=0; y(kp:kk)=0; 
e(1:kp-1)=0; e(kp:kk)=0; 
z(1:kp-1)=0; z(kp:kk)=0; % yzad
up(1:kk)=0;
ud(1:kk)=0;
ui(1:kk)=0;
u1(1:kk)=0;
u2(1:kk)=0;
z(10:kk) = 1;

T  = 2;
K  = 1.2;
Td = 5;
Ti = 20;
for k=kp:kk;
    %symulacja obiektu
    y(k)=0;
    for i=1:nb
        y(k)=y(k)+b(i)*u(k-i);
    end;
    for i=1:na
        y(k)=y(k)-a(i)*y(k-i);
    end;
    
    e(k)=z(k)-y(k);
    %e(k)=0.1;          % do wyznaczania czasu zdwojenia (PI)
    %e(k)=e(k-1)+0.1;   % do wyznaczania czasu wyprzedzenia (PD)
    
    % Podejscie pierwsze do wyznaczenia wartosci sterowania
    r2 = K*Td/T;                    
    r1 = K*(T/(2*Ti) - 2*Td/T - 1); 
    r0 = K*(1 + T/(2*Ti) + Td/T);   
    u1(k) = u(k-1)+r2*e(k-2)+r1*e(k-1)+r0*e(k);
    
    % Podejscie drugie do wyznaczenia wartosci sterowania
    up(k) = K*e(k);
    ud(k) = K*Td*(e(k)-e(k-1))/T;
    ui(k) = ui(k-1)+K/Ti*T*(e(k-1)+e(k))/2;
    u2(k) = up(k)+ud(k)+ui(k);
    
    u(k)=u1(k);
    %if(u(k)>umax); u(k) = umax; end;
    %if(u(k)<umin); u(k) = umin; end;
end;

%wyniki symulacji
figure; 
stairs((1:kk)*T,u,'b');
hold on;
stairs((1:kk)*T,up,'r');
stairs((1:kk)*T,ud+up,'k');
stairs((1:kk)*T,ui+ud+up,'g');
xlabel('k');
ylabel('u');

figure;
stairs((1:kk)*T,u1,'b');
xlabel('k');
ylabel('u');

figure; stairs(z,'k:');
hold on; stairs(y,'b');
xlabel('k');
ylabel('y_{zad}, y');
end

function [a,b,umin,umax,dumax,Tp] = obiekt1()
    Tp = 0.03;
    K=3; T=10;
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
    T2=0.1;
    [ld,md]=c2dm(K,[T1*T2 T1+T2 1],Tp,'zoh');
    b(1)=ld(2);
    b(2)=ld(3);
    a(1)=md(2);
    a(2)=md(3);
    umin=-1.0;
    umax=1.0;
    dumax=0.01;
end

function [a,b,umin,umax,dumax,Tp] = obiekt3()
    Tp = 0;
    a = [-1.6375, 0.67003];
    b = [0, 0, 0.035, 0.0307];
    umin = -Inf;
    umax = Inf;
    dumax = Inf;
end