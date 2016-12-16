cn = [0.75,0.00,0.00];
cm = [0.00,0.75,0.00];
cp = [0.00,0.00,0.75];
close all;
T  = 2;
K  = 1.2;
Ti = 20;
Td = 5;
[u,y,z]=PID_tests(T,K,Ti,Td,Inf);

[up1,yp1,zp1]=PID_tests(T,K-0.5,Ti,Td,Inf);
[up2,yp2,zp2]=PID_tests(T,K+0.5,Ti,Td,Inf);

[ui1,yi1,zi1]=PID_tests(T,K,Ti/2.5,Td,Inf);
[ui2,yi2,zi2]=PID_tests(T,K,Ti*10,Td,Inf);
% [ui1,yi1,zi1]=PID_tests(T,K,Ti*0.75,Td,Inf);
% [ui2,yi2,zi2]=PID_tests(T,K,Ti*2,Td,Inf);

[ud1,yd1,zd1]=PID_tests(T,K,Ti,Td-4,Inf);
[ud2,yd2,zd2]=PID_tests(T,K,Ti,Td+4,Inf);

component = 'd';
up = eval(sprintf('u%c2',component)); un = u; um = eval(sprintf('u%c1',component));
yp = eval(sprintf('y%c2',component)); yn = y; ym = eval(sprintf('y%c1',component));
if(component == 'p')
    p = 'K=1.7'; n = 'K=1.2'; m = 'K=0.7';
elseif(component == 'i')
    p = 'Ti=200'; n = 'Ti=20'; m = 'Ti=8';
    %p = 'Ti=40'; n = 'Ti=20'; m = 'Ti=8';
elseif(component == 'd')
    p = 'Td=9'; n = 'Td=5'; m = 'Td=1';
end

%wyniki symulacji
figure;
t = (1:length(u));
%t = 1:50;
hold on;
stairs(t,up(t),'Color', cp);
stairs(t,un(t),'Color', cn);
stairs(t,um(t),'Color', cm);
xlabel('k');
ylabel('u');
legend({p,n,m});

figure; 
hold on;
stairs(t,z  (t),'k:');
stairs(t,yp(t),'Color', cp);
stairs(t,yn(t),'Color', cn);
stairs(t,ym(t),'Color', cm);
xlabel('k');
ylabel('y_{zad}, y');
legend({'zadana',p,n,m});

close all;

figure;
t = (1:length(u));
hold on;
stairs(t,un(t),'b');
xlabel('k');
ylabel('u');

figure; 
hold on;
stairs(t,z  (t),'k:');
stairs(t,yn(t),'b');
xlabel('k');
ylabel('y_{zad}, y');