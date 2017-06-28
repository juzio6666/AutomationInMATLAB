global Hs Hp
delete(instrfindall); % zamkniecie wszystkich polaczen szeregowych
%clear all;
close all;
s = serial('COM4'); % COM9 to jest port utworzony przez mikrokontroler
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
while length(y)~=1100    % zbieramy 1000 pomiarow
    txt = fread(s,117); % odczytanie z portu szeregowego
                        % txt powinien zawieraæ "Y=%4d;U=%4d;"
                        % czyli np. "Y=1234;U=3232;"
    eval(char(txt'));   % wykonajmy to co otrzymalismy
    y=[y;Y];            % powiekszamy wektor y o element Y
    u=[u;U];            % powiekszamy wektor u o element U
    z=[z;Z];            % powiekszamy wektor z o element Z
end
dir = 'D:\Git\GitHub\LaTeX\CZASOPISMO2_Artyku³\czas_iteracji_vs_Hp_Hs\GPC2x2\';
name = sprintf('tp50ms_lambda1_Hp%d_Hs%d', Hp, Hs);

close all;
figure(1);   plot((0:(length(y)-1))*Tp,y(:,1),'r'); % wyswietlamy y w czasie
hold on;     plot((0:(length(y)-1))*Tp,y(:,2),'b'); % wyswietlamy y w czasie
hold on;   stairs((0:(length(y)-1))*Tp,z(:,1),'r--'); % wyswietlamy y w czasie
hold on;   stairs((0:(length(y)-1))*Tp,z(:,2),'b--'); % wyswietlamy y w czasie
ylabel('wartoœæ sygna³u wyjœciowego z zakresu [-2048, 2047]');
xlabel('czas (s)');
legend('wyjœcie 1', 'wyjœcie 2', 'zadana 1', 'zadana 2');
title(sprintf('Regulacja GPC2x2: Tp = %.2fs, Lambda = 1, Hp = %d, Hs = %d', Tp, Hp, Hs));
ylim([-2100,2100]);
savefig(1, [dir, 'fig\', name,'_output']);
saveas(1, [dir, 'png\', name,'_output'], 'png');

figure(2); stairs((0:(length(u)-1))*Tp,u(:,1),'r'); % wyswietlamy u w czasie
hold on;   stairs((0:(length(u)-1))*Tp,u(:,2),'b'); % wyswietlamy u w czasie
ylabel('wartoœæ sygna³u wejœciowego z zakresu [-2048, 2047]');
xlabel('czas (s)');
legend('sygna³ steruj¹cy 1', 'sygna³ steruj¹cy 2');
title(sprintf('Regulacja GPC2x2: Tp = %.2fs, Lambda = 1, Hp = %d, Hs = %d', Tp, Hp, Hs));
ylim([-2100,2100]);
savefig(2, [dir, 'fig\', name,'_input']);
saveas(2, [dir, 'png\', name,'_input'], 'png');


f = fopen([dir, 'csv\', name, '.txt'],'w');
fprintf(f,'k, t, u1, u2, y1, y2, z1, z2\n');
for i=1:length(u)
    fprintf(f,'%d, %.2f, %.2f, %.2f, %.0f, %.0f, %.0f, %.0f\n',i, (i-1)*Tp, u(i,1), u(i,2), y(i,1), y(i,2), z(i,1), z(i,2));
end
fclose(f);