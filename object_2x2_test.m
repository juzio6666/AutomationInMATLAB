clearvars;
%close all;

% Obiekt regulacji
%      input 1        input 2
Gs = [tf( 1,[.7 1]), tf( 5,[.3 1]);  % output 1
      tf( 1,[.5 1]), tf( 2,[.4 1])]; % output 2
Tp = 0.005;
Gz = c2d(Gs,Tp,'zoh');
ny = 2;
nu = 2;
% Based on STP script, pages 79 (tab. 3.1) and 86 (eq. 3.49)
for m=1:2
    for n=1:2
        alfa = exp(-Tp/Gs(m,n).Denominator{1}(1));
        B(m,n) = Gs(m,n).Numerator{1}(2)*(1-alfa);
        A(m,n) =                           -alfa ;
    end
end

% Y1/U1=G(1,1) and Y1/U2=G(1,2) => Y1 = G(1,1)*U1 + G(1,2)*U2
for m=1:2
    a(m,:)   = [A(m,1)+A(m,2), A(m,1)*A(m,2)];
    b(m,1,:) = [B(m,1), A(m,2)*B(m,1)];
    b(m,2,:) = [B(m,2), A(m,1)*B(m,2)];
end

kk = 1000;
y = zeros(2,kk);
u = [zeros(1,kk);ones(1,kk)];
for k = 1:kk
    for m=1:2
        if(k-2>=1)
            y(m,k) = b(m,1,1)*u(1,k-1) + b(m,1,2)*u(1,k-2) ...
                   + b(m,2,1)*u(2,k-1) + b(m,2,2)*u(2,k-2) ...
                   - a(m,1)  *y(m,k-1) - a(m,2)  *y(m,k-2);
        elseif(k-1>=1)
            y(m,k) = b(m,1,1)*u(1,k-1) + b(m,1,2)*0 ...
                   + b(m,2,1)*u(2,k-1) + b(m,2,2)*0 ...
                   - a(m,1)  *y(m,k-1) - a(m,2)  *0;
        else
            y(m,k) = b(m,1,1)*0        + b(m,1,2)*0 ...
                   + b(m,2,1)*0        + b(m,2,2)*0 ...
                   - a(m,1)  *0        - a(m,2)  *0;
        end            
    end
end

figure;
stairs((0:(kk-1))*Tp,y(1,:),'r');hold on;
stairs((0:(kk-1))*Tp,y(2,:),'b');hold off;