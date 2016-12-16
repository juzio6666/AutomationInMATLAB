clear all;
syms s z
licznik = expand((s+1)*(s+1));
mianownik = expand((s-2)*(s+3)*(s+4));

fprintf('Licznik: %s\n',char(licznik));
fprintf('Mianownik: %s\n',char(mianownik));

Gs = tf([1,2,1],[1,5,-2,-24]);

c2d(Gs, 0.5, 'zoh')
f = ((s+1)*(s+1)/((s-2)*(s+3)*(s+4)))/s;
sym A B C D
A = 3/20;
B = 
fprintf('%s\n',char(simplify(A/(s-2) + B/(s+3) + C/(s+4) + D/s)));

isequaln(A/(s-2) + B/(s+3) + C/(s+4) + D/s, 