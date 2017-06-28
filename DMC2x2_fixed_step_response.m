clear all;
load('DMC2x2_fixed_step_response');
y1_1 = y1(u1(:,1)~=0,:);
y1_1 = (y1_1(2:end,:)-ones(length(y1_1)-1,1)*y1_1(1,:))/1000;
y2_1 = y2(u2(:,2)~=0,:);
y2_1 = (y2_1(2:end,:)-ones(length(y2_1)-1,1)*y2_1(1,:))/500;
len = min(length(y1_1),length(y2_1));
len = 60;

s11 = y1_1(1:len,1);
s12 = y1_1(1:len,2);
s21 = y2_1(1:len,1);
s22 = y2_1(1:len,2);

for k = 1:len
    s(k,1,1) = s11(k);
    s(k,1,2) = s12(k);
    s(k,2,1) = s21(k);
    s(k,2,2) = s22(k);
end
figure;
hold on;
plot(s11);
plot(s12);
plot(s21);
plot(s22);
hold off;