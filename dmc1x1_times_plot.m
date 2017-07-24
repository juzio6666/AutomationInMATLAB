times = [  16   24   32   40   48   64   80  112  152  184  224  296  360  712 1060 1420];
D     = [  10   15   20   25   30   40   50   75  100  125  150  200  250  500  750 1000];

X=[ones(length(D),1),D'];
A=X\times';

plot(D,times); hold on;
x = D;
approx=[ones(size(x))',x']*A;
plot(x,approx); hold off;
csvwrite('../LaTeX/DUNNO_Pomiar_czasu_algorytmow_regulacji/dane/dmc1x1time.csv',[D' times' approx]);
