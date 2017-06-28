%% Skrypt do sprawdzania projektu 1 z STP

%% Treœæ zadania
% _Obiekt dynamiczny opisany jest transmitancj¹ ci¹g³¹_
%
% $$ G(s) = \frac{(s-z_0)(s-z_1)}{(s-b_0)(s-b_1)(s-b_2)} $$
%
% Treœci zadañ projektowych ró¿ni¹ siê zerami ($z_0$, $z_1$), biegunami 
% ($b_0$, $b_1$, $b_2$) i warunkiem pocz¹tkowym symulacji, np.:
zerac    = [-1, -10];      % [z_0, z_1]
biegunyc = [11, -12, -13]; % [b_0, b_1, b_2]
x0       = [-1 -2 5]';     % 

%% Zadanie 1
% _Wyznaczyæ transmitancjê dyskretn¹ $G(z)$ . Zastosowaæ okres próbkowania 
% $0.5$ sek. i ekstrapolator zerowego rzêdu. Okreœliæ zera i bieguny 
% transmitancji ci¹g³ej i dyskretnej._
%
% Rozwi¹zanie dotycz¹ce czêœci ci¹g³ej mo¿na odczytaæ z transmitancji z
% opisu zadania. Aby uczyniæ to samo dla transmitancji dyskretnej mo¿na
% albo przekszta³ciæ bieguny przy u¿yciu funkcji exponent, albo
% przekszta³ciæ transmitancjê ci¹g³¹ na dyskretn¹ i odczytaæ je z wyniku
Tp = 0.5;
Gs = tf(poly(zerac), poly(biegunyc));

Gz = c2d(Gs, Tp, 'zoh');

Kc = prod(-zerac)/prod(-biegunyc);
zerad = exp(zerac*Tp);
biegunyd = exp(biegunyc*Tp);
Kd = prod(1-zerad)/prod(1-biegunyd);
K  = Kc/Kd;
Gz2 = K*tf(poly(zerad),poly(biegunyd),Tp);

%% Zadanie 2

%% Zadanie 3

%% Zadanie 4

%% Zadanie 5

%% Zadanie 6
