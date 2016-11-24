clear all;
red   = [0.75, 0.00, 0.00];
green = [0.00, 0.75, 0.00];
blue  = [0.00, 0.00, 0.75];
T  = 2;
K  = 1.2;
Ti = 20;
Td = 5;
[u0,y0,z0,~]=PID_tests(T,K,Ti,Td,Inf);

dK = 0.5;
[u0,y0,z0,~]=PID_tests(T,K,Inf,0,Inf);
[u1,y1,~,~]=PID_tests(T,K-dK,Inf,0,Inf);
[u2,y2,~,~]=PID_tests(T,K+dK,Inf,0,Inf);

figure;
t = (1:length(u0));
hold on;
stairs(t,u2(t),'Color', blue);
stairs(t,u0(t),'Color', red);
stairs(t,u1(t),'Color', green);
xlabel('k');
ylabel('u');
legend(sprintf('K=%3.1f',K-dK),...
       sprintf('K=%3.1f',K),...       'spodziewane u(k) (z anti-windup)',...
       sprintf('K=%3.1f',K+dK)'...       'faktyczne u(k) (z anti-windup)'
);
   
figure; 
hold on;
stairs(t,z0(t),'k:');
stairs(t,y2(t),'Color',blue);
stairs(t,y0(t),'Color',red);
stairs(t,y1(t),'Color',green);
xlabel('k');
ylabel('y_{zad}, y');
legend('zadana',...
       sprintf('K=%3.1f',K-dK),...
       sprintf('K=%3.1f',K),...       'spodziewane u(k) (z anti-windup)',...
       sprintf('K=%3.1f',K+dK)'...       'faktyczne u(k) (z anti-windup)'
);


[u0,y0,z0,~]=PID_tests(T,K,Ti,0,Inf);
[u1,y1,~,~]=PID_tests(T,K,Ti*0.75,0,Inf);
[u2,y2,~,~]=PID_tests(T,K,Ti*2,0,Inf);

figure;
t = (1:length(u0));
hold on;
stairs(t,u2(t),'Color', blue);
stairs(t,u0(t),'Color', red);
stairs(t,u1(t),'Color', green);
xlabel('k');
ylabel('u');
legend(sprintf('Ti=%3.1f',Ti*2),...
       sprintf('Ti=%3.1f',Ti),...       'spodziewane u(k) (z anti-windup)',...
       sprintf('Ti=%3.1f',Ti*0.75)'...       'faktyczne u(k) (z anti-windup)'
);
   
figure; 
hold on;
stairs(t,z0(t),'k:');
stairs(t,y2(t),'Color',blue);
stairs(t,y0(t),'Color',red);
stairs(t,y1(t),'Color',green);
xlabel('k');
ylabel('y_{zad}, y');
legend('zadana',...
       sprintf('Ti=%3.1f',Ti*2),...
       sprintf('Ti=%3.1f',Ti),...       'spodziewane u(k) (z anti-windup)',...
       sprintf('Ti=%3.1f',Ti*0.75)'...       'faktyczne u(k) (z anti-windup)'
);




close all;

[u0,y0,z0,~]=PID_tests(T,K,Inf,Td,Inf);
[u1,y1,~,~]=PID_tests(T,K,Inf,Td-4,Inf);
[u2,y2,~,~]=PID_tests(T,K,Inf,Td+5,Inf);

figure;
t = (1:length(u0));
hold on;
stairs(t,u2(t),'Color', blue);
stairs(t,u0(t),'Color', red);
stairs(t,u1(t),'Color', green);
xlabel('k');
ylabel('u');
legend(sprintf('Td=%3.1f',Td-4),...
       sprintf('Td=%3.1f',Td),...       'spodziewane u(k) (z anti-windup)',...
       sprintf('Td=%3.1f',Td+5)'...       'faktyczne u(k) (z anti-windup)'
);
   
figure; 
hold on;
stairs(t,z0(t),'k:');
stairs(t,y2(t),'Color',blue);
stairs(t,y0(t),'Color',red);
stairs(t,y1(t),'Color',green);
xlabel('k');
ylabel('y_{zad}, y');
legend('zadana',...
       sprintf('Td=%3.1f',Td-4),...
       sprintf('Td=%3.1f',Td),...       'spodziewane u(k) (z anti-windup)',...
       sprintf('Td=%3.1f',Td+5)'...       'faktyczne u(k) (z anti-windup)'
);