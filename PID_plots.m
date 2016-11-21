cn = [0.75,0.00,0.00];
cm = [0.00,0.75,0.00];
cp = [0.00,0.00,0.75];

T  = 2;
K  = 1.2;
Ti = 20;
Td = 5;
[u,y,z]=PID_tests(T,K,Ti,Td);

[up1,yp1,zp1]=PID_tests(T,K-0.5,Ti,Td);
[up2,yp2,zp2]=PID_tests(T,K+0.5,Ti,Td);

[ui1,yi1,zi1]=PID_tests(T,K,Ti/2.5,Td);
[ui2,yi2,zi2]=PID_tests(T,K,Ti*10,Td);

[ud1,yd1,zd1]=PID_tests(T,K,Ti,Td-4);
[ud2,yd2,zd2]=PID_tests(T,K,Ti,Td+4);

component = 'i';
up = eval(sprintf('u%c2',component)); un = u; um = eval(sprintf('u%c1',component));
yp = eval(sprintf('y%c2',component)); yn = y; ym = eval(sprintf('y%c1',component));
if(component == 'p')
    p = 'K=1.7'; n = 'K=1.2'; m = 'K=0.7';
elseif(component == 'i')
    p = 'Ti=8'; n = 'Ti=20'; m = 'Ti=200';
elseif(component == 'd')
    p = 'Td=1'; n = 'Td=5'; m = 'Td=9';
end

%wyniki symulacji
figure;
t = (1:length(u1));
%t = 1:50;
hold on;
stairs(t*T,up(t),'Color', cp);
stairs(t*T,un(t),'Color', cn);
stairs(t*T,um(t),'Color', cm);
xlabel('k');
ylabel('u');
legend({p,n,m});

figure; 
hold on;
stairs(t*T,z  (t),'k:');
stairs(t*T,yp(t),'Color', cp);
stairs(t*T,yn(t),'Color', cn);
stairs(t*T,ym(t),'Color', cm);
xlabel('k');
ylabel('y_{zad}, y');
legend({'zadana',p,n,m});

close all;

figure;
t = (1:length(u1));
hold on;
stairs(t*T,un(t),'b');
xlabel('k');
ylabel('u');

figure; 
hold on;
stairs(t*T,z  (t),'k:');
stairs(t*T,yn(t),'b');
xlabel('k');
ylabel('y_{zad}, y');