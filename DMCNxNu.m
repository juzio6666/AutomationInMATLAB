function DMCNxNu(S,N,Nu)
    object = @object2;
    
    kmax = 1000;
    kstart = 10;
    ny = 2;
    nu = 2;

    D  = size(S,3); 
    
    yzad = [10, 5]'*ones(1,kmax+N+1);

    nY = N*ny;
    ndU = Nu*nu;
    y = zeros(ny,kmax);
    du = zeros(nu,kmax);
    u = zeros(nu,kmax);
        
    M = zeros(nY, ndU);
    for i = 1:N
        for j=1:min(i,Nu)
            M(((i-1)*ny+1):(i*ny),(((j-1)*nu+1):(j*nu))) = S(i-j+1,:,:);
        end
    end
    
    Mp = zeros(ny*N, nu*(D-1));
    for i = 1:N
        for j = 1:(D-1)
            tmp = S(D,:,:);
            if(i+j<D)
                tmp = S(i+j,:,:);
            end
            Mp(((i-1)*ny+1):(i*ny),((j-1)*nu+1):(j*nu)) = tmp-S(j,:,:);
        end
    end
    
    lambda = 100;
    Psi_ = diag(ones(ny*N,1));
    Lambda_ = diag(lambda*ones(nu*Nu,1));
    
    K = (M'*Psi_*M+Lambda_)^(-1)*M'*Psi_;
    dUp = zeros(nu*(D-1),1);
    Y   = zeros(ny*N,1);
    Yzad= zeros(ny*N,1);
    
    
    
    for k = kstart:kmax
        if(k-1 > 0)
            y = object(u(:,:), y(:,:),k);
        else
            y(:,k) = zeros(1,ny);
        end
        if k < kstart
            du(:,k) = 0;
            if(k-1 > 0)
                u(:,k) = u(:,k-1)+du(:,k);
            else
                u(:,k) = zeros(1,nu);
            end
            continue
        end
        
        for i = 1:(D-1)
            dUp(((i-1)*nu+1):(i*nu)) = du(:,(k-i));
        end
        for j = 1:N
            Yzad(((j-1)*ny+1):(j*ny)) = yzad(:,(k+j));
            Y(((j-1)*ny+1):(j*ny)) = y(:,k)';
        end
        
        Y0 = Y+Mp*dUp;
        dU = K*(Yzad - Y0); % <-- wyznaczamy nowe wartoœci sygna³u steruj¹cego
        du(:,k) = dU(1:nu);
        u(:,k) = u(:,k-1)+du(:,k);
    end
    plot(1:kmax, y);
end

function [y] = object1(u,y)
    y(1) = y(1)*0.9+u(1)*0.1;
    y(2) = y(2)*0.8+u(2)*0.2;
end

function [y] = object2(u,y,k)
    ilwe=2; ilwy=2;
    Tp=0.03;%okres próbkowania
    K=1; T=0.7; alfa=exp(-Tp/T); B11=K*(1-alfa); A11=-alfa; %G11=1/0.7s+1
    K=5; T=0.3; alfa=exp(-Tp/T); B12=K*(1-alfa); A12=-alfa; %G12=5/0.3s+1    
    K=1; T=0.5; alfa=exp(-Tp/T); B21=K*(1-alfa); A21=-alfa; %G21=1/0.5s+1    
    K=2; T=0.4; alfa=exp(-Tp/T); B22=K*(1-alfa); A22=-alfa; %G22=2/0.4s+1
    
    %sprowadzenie do wspólnego mianownika (równania ró¿nicowe)
    na=2; nb=2;
    %a(m,i), m - wyjœcie, i - przesuniêcie
    a(1,1)=A11+A12; a(1,2)=A11*A12;
    a(2,1)=A21+A22; a(2,2)=A21*A22;
    %b(m,n,i), m - wyjœcie, n - wejœcie, i - przesuniêcie
    b(1,1,1)=B11; b(1,1,2)=B11*A12;
    b(1,2,1)=B12; b(1,2,2)=B12*A11;
    b(2,1,1)=B21; b(2,1,2)=B21*A22;
    b(2,2,1)=B22; b(2,2,2)=B22*A21;
        k
    for m=1:ilwy;
        y(m,k)=0.0;
        for n=1:ilwe;
            for i=1:nb;
                y(m,k)=y(m,k)+b(m,n,i)*u(n,k-i);
            end;
        end;
        for i=1:na;
            y(m,k)=y(m,k)-a(m,i)*y(m,k-i);
        end;
    end;
end