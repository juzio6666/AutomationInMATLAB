T  = 0.03;
K  = 1.0;
Ti = 0.1;
Td = 0.0;
umin = -1.0;
umax =  1.0;
[u0,y0,z0,awup0]=PID_tests(T,K,Ti,Td,Inf);
[u1,y1,z1,awup1]=PID_tests(T,K,Ti,Td,0.1);

figure;
t = (1:length(u));
hold on;
stairs(t*T,ones(1,length(t))*umax,'k-.');
stairs(t*T,u0(t),'b--');
stairs(t*T,u1(t),'r--');
stairs(t*T,max(min(umax,u0(t)),umin),'b');
stairs(t*T,max(min(umax,u1(t)),umin),'r');
xlabel('k');
ylabel('u');
legend('ograniczenie umax=1.0',...
       'spodziewane u(k) (bez anti-windup)',...
       'spodziewane u(k) (z anti-windup)',...
       'faktyczne u(k) (bez anti-windup)',...
       'faktyczne u(k) (z anti-windup)');

figure; 
hold on;
stairs(t*T,z1(t),'k:');
stairs(t*T,y0(t),'b');
stairs(t*T,y1(t),'r');
xlabel('k');
ylabel('y_{zad}, y');
legend('zadana',...
       'y(k) (bez anti-windup)',...
       'y(k) (z anti-windup)');