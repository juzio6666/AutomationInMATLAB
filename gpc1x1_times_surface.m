% na =      1    2    3    4    5
times = [ 5.2  6.0  6.6  7.3  8.0; % nb = 1
          5.7  6.5  7.3  7.9  8.7; % nb = 2
          6.5  7.2  7.9  8.6  9.2; % nb = 3
          7.1  7.8  8.5  9.1  9.8; % nb = 4
          7.7  8.4  9.1  9.8 10.5; % nb = 5
          8.4  9.0  9.7 10.3 11.1; % nb = 6
          9.0  9.6 10.3 11.0 11.6; % nb = 7
%        nan  nan  nan  nan  nan ; % nb = 8
%        nan  nan  nan  nan  nan ; % nb = 9
         10.7 11.4 12.1 12.8 13.5];% nb = 10
     
na = 1:5;
nb = [1:7,10];

[NA,NB]=meshgrid(na,nb);
surface(NA,NB,times);