function DMCNxNu(S,N,Nu)
    u1 = [zeros(2,100), [ones(1,900); zeros(1,900)]];
    u2 = [zeros(2,100), [zeros(1,900); ones(1,900)]];
    y1 = 0*u1;
    y2 = 0*u2;
    for k = 10:1000
        y1(:,k) = object(u1(:,k-1), y1(:,k-1));
        y2(:,k) = object(u2(:,k-1), y2(:,k-1));
    end
    
    kmax = 1000;
    kstart = 200;
    N = 4;
    Nu = 2;
    S(:,1,1) = y1(1,100:200);
    S(:,1,2) = y1(2,100:200);
    S(:,2,1) = y2(1,100:200);
    S(:,2,2) = y2(2,100:200);
    yzad = [10, 5]'*ones(1,kmax+N+1);
    
    D = size(S,1);
    ny = size(S,2);
    nu = size(S,3);
    for i = (D+1):(D-1+N)
        S(i,:,:) = S(D,:,:);
    end
    nY = N*ny;
    ndU = Nu*nu;
    y = zeros(ny,0);
    du = zeros(nu,0);
    u = zeros(nu,0);
        
    M = zeros(nY, ndU);
    for i = 1:N
        for j=1:min(i,Nu)
            M(((i-1)*ny+1):(i*ny),(((j-1)*nu+1):(j*nu))) = S(i-j+1,:,:);
        end
    end
    
    Mp = zeros(ny*N, nu*(D-1));
    for i = 1:N
        for j = 1:(D-1)
            Mp(((i-1)*ny+1):(i*ny),((j-1)*nu+1):(j*nu)) = S(i+j,:,:)-S(j,:,:);
        end
    end
    
    lambda = 12;
    Psi_ = diag(ones(ny*N,1));
    Lambda_ = diag(lambda*ones(nu*Nu,1));
    
    K = (M'*Psi_*M+Lambda_)^(-1)*M'*Psi_;
    dUp = zeros(nu*(D-1),1);
    Y   = zeros(ny*N,1);
    Yzad= zeros(ny*N,1);
    
    for k = 1:kmax
        if(k-1 > 0)
            y(:,k) = object(u(:,k-1), y(:,k-1));
        else
            y(:,k) = [1,2];
        end
        if k < kstart
            du(:,k) = 0;
            if(k-1 > 0)
                u(:,k) = u(:,k-1)+du(:,k);
            else
                u(:,k) = [0.1,0.2];
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

function [y] = object(u,y)
    y(1) = y(1)*0.9+u(1)*0.1;
    y(2) = y(2)*0.8+u(2)*0.2;
end