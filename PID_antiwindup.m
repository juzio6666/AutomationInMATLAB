T  = 0.03;
K  = 1.0;
Ti = 0.1;
Td = 2;
[u1,y1,z1,awup1]=PID_tests(T,K,Ti,Td,0.01);

figure;
t = (1:length(u));
hold on;
stairs(t*T,u1(t),'b');
xlabel('k');
ylabel('u');

figure; 
hold on;
stairs(t*T,z1(t),'k:');
stairs(t*T,y1(t),'b');
xlabel('k');
ylabel('y_{zad}, y');