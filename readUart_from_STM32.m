delete(instrfindall); % zamkniecie wszystkich polaczen szeregowych
clearvars;

s = serial('COM3'); % COM9 to jest port utworzony przez mikrokontroler
set(s,'BaudRate',115200);
set(s,'StopBits',1);
set(s,'Parity','none');
set(s,'DataBits',8); 
set(s,'Timeout',100000);
set(s,'InputBufferSize',1000);
set(s,'Terminator',13);
fopen(s); % otwarcie kanalu komunikacyjnego

Tp = 0.005; % czas z jakim probkuje regulator
y = [];    % wektor wyjsc obiektu
z = [];    % wektor wyjsc obiektu6
u = [];    % wektor wejsc (sterowan) obiektu

while length(y)~=2000    % zbieramy 100 pomiarow
    txt = fread(s,84);  % odczytanie z portu szeregowego
                        % txt powinien zawieraæ "Y=%4d;U=%4d;"
                        % czyli np. "Y=1234;U=3232;"
    %disp(char(txt'));disp('\n');
    eval(char(txt'));   % wykonajmy to co otrzymalismy
    y=[y;Y];            % powiekszamy wektor y o element Y
    z=[z;Z];            % powiekszamy wektor y o element Y
    u=[u;U];            % powiekszamy wektor u o element U
end

%close all;
figure;
subplot(2,1,1);
plot(y); hold on;
plot(z,'--'); hold off;
subplot(2,1,2);
stairs(u); 