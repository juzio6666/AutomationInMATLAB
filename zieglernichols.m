clear all;
close all;
T = 0.01;
Wzm=20.0; 
T1=1; 
T2=0.3;
[ld,md]=c2dm(Wzm,[T1*T2 T1+T2 1],T,'zoh');
b(1)=ld(2);
b(2)=ld(3);
a(1)=md(2);
a(2)=md(3);
na=length(a); nb=length(b); 

kp=5; kk=1000;

u(1:kk)=0;
y(1:kk)=0;
e(1:kk)=0;
z(1:kk)=0; % yzad
up(1:kk)=0;
ud(1:kk)=0;
ui(1:kk)=0;

z(kp:kk) = 250;

%% parametry do wprowadzenia w oscylacje
Ku = 13.1; % wzmocnienie krytyczne
Tu = 0.215; % okres oscylacji biorac pod uwage czas probkowania 0.01s

%% parametry PID wyznaczone na podstawie Z-N z polskiej wikipedii
typ_reg = 'PID';
if strcmp(typ_reg, 'P')
    K = 0.5*Ku;
    Ti = Inf;
    Td = 0;
elseif strcmp(typ_reg, 'PI')
%     K = 0.45*Ku;
%     Ti = Tu/1.2;
    K = Ku/3.2;
    Ti = 2.2*Tu;
    Td = 0;
elseif strcmp(typ_reg, 'PID')
%     K = 0.6*Ku;
%     Ti = 0.5*Tu;
%     Td = 0.125*Tu;
    K = Ku/2.2;
    Ti = 2.2*Tu;
    Td = Tu/6.3;
else
    K = 13.1;
    Ti = Inf;
    Td = 0;
end

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
    
    up(k) = K*e(k);
    ud(k) = K*Td*(e(k)-e(k-1))/T;
    ui(k) = ui(k-1)+K/Ti*T*(e(k-1)+e(k))/2;
    u(k)  = up(k)+ud(k)+ui(k);
end;

figure; plot((0:length(y)-1)*T, y);
figure; plot((0:length(u)-1)*T, u);