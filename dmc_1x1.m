function [du, K, Ku, Ke, Mp] = dmc_1x1(S,N,Nu,Lambda,Psi,dUp,Y,Yzad)
[K, Mp, Ku, Ke] = dmc_1x1_init(S,N,Nu,Lambda,Psi);
%du = dmc_1x1_loop(N,length(S),K,Mp,dUp,Y,Yzad);
du = dmc_1x1_loop_min(length(S),Ku,Ke,dUp,Y(1),Yzad(1));
end

function [K, Mp, Ku, Ke] = dmc_1x1_init(S,N,Nu,Lambda,Psi)
    D = length(S);
    %% Matrices, vectors and scalars preparation
    % Memmory assignment
    M = zeros(N,Nu);
    Mp = zeros(N,D-1);

    % Matrix M
    for row = 1:N
       for col = 1:Nu
            if(row-col+1 >= 1)
                M(row,col) = S(row-col+1);
            end
       end
    end

    % Matrix Mp
    for row = 1:N
       for col = 1:(D-1)
            Mp(row,col) = S(min(row+col,D)) - S(col);
       end
    end
    K = (M'*Psi*M+Lambda)^(-1)*M';
    Ke = sum(K(1,:));
    Ku = K(1,:)*Mp;    
end

%% Algorithm itself (extensive form)
function du = dmc_1x1_loop(N,D,K,Mp,dUp,Y,Yzad)
    Y0 = Y+Mp*dUp;
    dU = K*(Yzad-Y0);
    du = dU(1);
end

%% Algorithm itself (minimalistic form assuming Yzad constant on the prediction horison)
function du = dmc_1x1_loop_min(D,Ku,Ke,dUp,y,yzad)
    du = Ke*(yzad-y) - Ku*dUp;
end
