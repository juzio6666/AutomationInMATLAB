delete(instrfindall); % zamkniecie wszystkich polaczen szeregowych
%clear all;
close all;
s = serial('COM3'); % COM9 to jest port utworzony przez mikrokontroler
set(s,'BaudRate',115200);
set(s,'StopBits',1);
set(s,'Parity','none');
set(s,'DataBits',8);
set(s,'Timeout',1);
set(s,'InputBufferSize',1000);
set(s,'Terminator',13);
fopen(s); % otwarcie kanalu komunikacyjnego

Tp = 0.05;
y = [];    % wektor wyjsc obiektu
u = [];    % wektor wejsc (sterowan) obiektu
z = [];    % wektor wartosci zadanych
while length(y)~=3000    % zbieramy 1000 pomiarow
    txt = fread(s,62+35); % odczytanie z portu szeregowego
                        % txt powinien zawieraæ "Y=%4d;U=%4d;"
                        % czyli np. "Y=1234;U=3232;"
    disp(char(txt'));
    eval(char(txt'));   % wykonajmy to co otrzymalismy
    y=[y;[y1 y2]];            % powiekszamy wektor y o element Y
    u=[u;[u1 u2]];            % powiekszamy wektor u o element U
    z=[z;[z1 z2]];            % powiekszamy wektor u o element U
end
close all;
figure(1); 
subplot(2,1,1);
plot((0:(length(y)-1))*Tp,y); % wyswietlamy y w czasie
hold on; 
stairs((0:(length(z)-1))*Tp,z); % wyswietlamy u w czasie
legend('y_1', 'y_2', 'y_1^{zad}', 'y_2^{zad}');
ylabel('y, y^{zad}');
xlabel('czas (s)');
subplot(2,1,2);
stairs((0:(length(u)-1))*Tp,u); % wyswietlamy u w czasie
ylabel('u');
legend('u_1', 'u_2');
xlabel('czas (s)');