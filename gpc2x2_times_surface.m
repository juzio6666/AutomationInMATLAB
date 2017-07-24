% na =      1    2    3    4    5
times = [33.6 36.8 41.6 46.4 51.2; % nb = 1
         38.4 43.2 46.4 51.2 56.0; % nb = 2
         43.2 48.0 51.2 56.0 60.8; % nb = 3
         46.4 51.2 56.0 59.2 64.0; % nb = 4
         51.2 56.0 60.8 65.6 70.4; % nb = 5
         56.0 60.8 64.0 70.4 73.6; % nb = 6
         60.8 64.0 68.8 73.6 78.4; % nb = 7
%        nan  nan  nan  nan  nan ; % nb = 8
%        nan  nan  nan  nan  nan ; % nb = 9
         73.6 78.4 81.6 88.0 91.2];% nb = 10
     
na = 1:5;
nb = [1:7,10];

[NA,NB]=meshgrid(na,nb);
surface(NA,NB,times);
csvwrite('../LaTeX/DUNNO_Pomiar_czasu_algorytmow_regulacji/dane/gpc2x2time.csv',[NA, NB, times]);