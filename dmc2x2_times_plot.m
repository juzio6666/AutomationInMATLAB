times = [ 112  168  216  272  320  424  528  800 1050 1300 1580 1980 2470];
D     = [  10   15   20   25   30   40   50   75  100  125  150  200  250];

X=[ones(length(D),1),D'];
A=X\times';

plot(D,times); hold on;
x = D;
approx=[ones(size(x))',x']*A;
plot(x,approx); hold off;
csvwrite('../LaTeX/DUNNO_Pomiar_czasu_algorytmow_regulacji/dane/dmc2x2time.csv',[D' times' approx]);
