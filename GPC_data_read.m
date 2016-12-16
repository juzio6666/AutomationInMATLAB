global Hs Hp
delete(instrfindall); % zamkniecie wszystkich polaczen szeregowych
%clear all;
close all;
s = serial('COM11'); % COM9 to jest port utworzony przez mikrokontroler
set(s,'BaudRate',115200);
set(s,'StopBits',1);
set(s,'Parity','none');
set(s,'DataBits',8);
set(s,'Timeout',100);
set(s,'InputBufferSize',1000);
set(s,'Terminator',13);
fopen(s); % otwarcie kanalu komunikacyjnego

Tp = 0.05; % czas z jakim probkuje regulator
e = [];    % wektor uchybow
y = [];    % wektor wyjsc obiektu
u = [];    % wektor wejsc (sterowan) obiektu
z = [];
while length(y)~=500    % zbieramy 1000 pomiarow
    txt = fread(s,117); % odczytanie z portu szeregowego
                        % txt powinien zawieraæ "Y=%4d;U=%4d;"
                        % czyli np. "Y=1234;U=3232;"
    eval(char(txt'));   % wykonajmy to co otrzymalismy
    y=[y;Y];            % powiekszamy wektor y o element Y
    u=[u;U];            % powiekszamy wektor u o element U
    z=[z;Z];            % powiekszamy wektor z o element Z
end
close all;
figure(1); plot((0:(length(y)-1))*Tp,y); % wyswietlamy y w czasie
hold on; stairs((0:(length(u)-1))*Tp,u); % wyswietlamy u w czasie
ylabel('wartoœæ sygna³u wejœciowego/wyjœciowego z zakresu [-2048, 2047]');
xlabel('czas (s)');
legend('wyjœcie obiektu', 'sygna³ steruj¹cy');
fprintf('Hp: %d, Hs: %d\n', Hp, Hs);
czas = input('Jaki by³ czas? ','s');
title(sprintf('Regulacja GPC2x2: Tp = %.2fs, Lambda = 15, Hp = %d, Hs = %d (czas iteracji = %s)', Tp, Hp, Hs, czas));
ylim([-2000,2000]);
dir = 'D:\Git\GitHub\LaTeX\CZASOPISMO2_Artyku³\czas_iteracji_vs_Hp_Hs\GPC2x2\';
name = sprintf('tp50ms_lambda15_Hp%d_Hs%d', Hp, Hs);
savefig([dir, 'fig\', name]);
saveas(1, [dir, 'png\', name], 'png');
f = fopen([dir, 'csv\', name, '.txt'],'w');
fprintf(f,'k, t, u, y, z\n');
for i=1:length(u)
    fprintf(f,'%d, %.2f, %.2f, %.0f, %.0f\n',i, (i-1)*Tp, u(i), y(i), z(i));
end
fclose(f);
close all;