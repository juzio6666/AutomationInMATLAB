global Hs Hp
delete(instrfindall); % zamkniecie wszystkich polaczen szeregowych
%clear all;
close all;
s = serial('COM16'); % COM9 to jest port utworzony przez mikrokontroler
set(s,'BaudRate',115200);
set(s,'StopBits',1);
set(s,'Parity','none');
set(s,'DataBits',8);
set(s,'Timeout',1);
set(s,'InputBufferSize',1000);
set(s,'Terminator',13);
fopen(s); % otwarcie kanalu komunikacyjnego

Tp = 0.05; % czas z jakim probkuje regulator
y = [];    % wektor wyjsc obiektu
u = [];    % wektor wejsc (sterowan) obiektu
odp = [];
k = 1;
odpk1 = Inf;
odpk2 = Inf;
D = 100;
while k<800   % zbieramy 1000 pomiarow
    txt = fread(s,50);  % odczytanie z portu szeregowego
                        % txt powinien zawieraæ "Y=%4d;U=%4d;"
                        % czyli np. "Y=1234;U=3232;"                        
    eval(char(txt'));   % wykonajmy to co otrzymalismy
    %char(txt')
    y(k,1) = y1;
    y(k,2) = y2;
    u(k,1) = u1;
    u(k,2) = u2;
    if(odpk1 == Inf && u(k,1) > 5)
       odpk1 = k+1;
    end
    if(odpk2 == Inf && u(k,2) > 5)
       odpk2 = k+1;
    end
    fprintf('odpk1 = %d; odpk2 = %d;\n', odpk1, odpk2);
end
close all;
figure(1); plot((0:(length(y)-1))*Tp,y); % wyswietlamy y w czasie
hold on; stairs((0:(length(u)-1))*Tp,u); % wyswietlamy u w czasie

figure(2);
plot((0:(D-1))*Tp,(y(odpk1+(0:(D-1)),1)-y(odpk1-1,1))/1000);
hold on;
plot((0:(D-1))*Tp,(y(odpk1+(0:(D-1)),2)-y(odpk1-1,2))/1000);
plot((0:(D-1))*Tp,(y(odpk2+(0:(D-1)),1)-y(odpk2-1,1))/200);
plot((0:(D-1))*Tp,(y(odpk2+(0:(D-1)),2)-y(odpk2-1,2))/200);