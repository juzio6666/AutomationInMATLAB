function [du, Knu] = gpc_1x1(a,b,N,Nu,Lambda,Psi,u,y,Yzad)
%%
% y = y(k-na), y(k-na+1), ... , y(k-2), y(k-1), y(k) -- d³ugoœæ na+1
% u = u(k-nb), u(k-na+1), ... , u(k-2), u(k-1),      -- d³ugoœæ nb
%
% y(k) --> y(na+1)
% u(k-1) --> u(nb)

Knu = gpc_1x1_const_data(a,b,N,Nu,Lambda,Psi);
du = gpc_1x1_loop(a,b,N,Knu,u,y,Yzad);

end

function Knu = gpc_1x1_const_data(a,b,N,Nu,Lambda,Psi)
    %% Here should be all of the correctness checking.
    % N check
    if(~isscalar(N)); error('N should ba a scalar'); end

    % Nu check
    if(~isscalar(Nu)); error('Nu should ba a scalar'); end

    % N >= Nu check
    if(N < Nu)
        warning('N should not be lower than Nu -- setting Nu to N');
        Nu=N;
    end

    % Lambda check
    if(isvector(Lambda))
        if(length(Lambda) ~= Nu)
            error('Lambda vector is of length %d instead of Nu=%d',length(Lambda),Nu);
        else % OK -- vector of diagonal elements
            Lambda = diag(Lambda);  
        end
    elseif(~isdiag(Lambda)); error('Lambda should be a diagonal matrix');
    else
        if(size(Lambda,1)~= Nu); error('Lambda should be a diagonal matrix with size NuxNu=%dx%d instead of %dx%d',Nu,Nu,size(Lambda));
        else % OK -- diagonal matrix
        end
    end
    if(min(diag(Lambda)) < 0 ); error('Diagonal of Lambda should contain only positive values'); end

    % Psi check 
    if(isvector(Psi))
        if(length(Psi) ~= N)
            error('Psi vector is of length %d instead of N=%d',length(Psi),N);
        else % OK -- vector of diagonal elements
            Psi = diag(Psi);
        end
    elseif(~isdiag(Psi)); error('Psi should be a diagonal matrix');
    else 
        if(size(Psi,1)~= N); error('Psi should be a diagonal matrix with size NxN=%dx%d instead of %dx%d',N,N,size(Psi));
        else % OK -- diagonal matrix
        end
    end
    if(min(diag(Psi)) <= 0 ); error('Diagonal of Psi should contain only nonnegative values'); end

    %% Matrices, vectors and scalars preparation
    % Memmory assignment
    M = zeros(N,Nu);
    S = zeros(N,1);
    
    nb = length(b);
    na = length(a);
    
    % Vector S
    for j=1:N
        S(j) = 0;
        for i=1:min(j,nb); S(j) = S(j) + b(i); end
        for i=1:min(j-1,na); S(j) = S(j) - a(i)*S(j-i); end
    end
    
    % Matrix M
    for row = 1:N
       for col = 1:Nu
            if(row-col+1 >= 1)
                M(row,col) = S(row-col+1);
            end
       end
    end

    K = (M'*Psi*M+Lambda)^(-1)*M';    
    Knu = K(1,:);
end

%% Algorithm itself (extensive form)
function du = gpc_1x1_loop(a,b,N,Knu,u,y,Yzad)
    % Yzad check
    if(isscalar(Yzad)); Yzad=ones(N,1)*Yzad;
    elseif(isvector(Yzad))
        if(length(Yzad) ~= N); error('Wrong length of Yzad vector -- got %d expected N=%d',length(Yzad),N);
        elseif(size(Yzad,1) == 1); Yzad=Yzad'; 
        else % Yzad is OK as is
        end
    else; error('Yzad should be a 1xN = 1x%d vector',N); 
    end
    
    nb = length(b);
    na = length(a);
    
    % u check
    if(isvector(u))
        if(length(u) ~= nb); error('Vector u should have a size of nb=%d instead of %d',nb,length(u));
        else % u is OK as is
        end
    else
        error('u is not a vector');
    end    
    
    % y check
    if(isvector(y))
        if(length(y) ~= (na+1)); error('Vector y should have a size of na+1=%d instead of %d',(na+1),length(y));
        else % u is OK as is
        end
    else
        error('y is not a vector');
    end    
    
    ym = 0;
    for i=1:nb; ym = ym + b(i)*u(nb+1-i); end % u(k-i) -k=nb+1-> u(nb+1-i)
    for i=1:na; ym = ym - a(i)*y(na+1-i); end % y(k-i) -k=na+1-> y(na+1-i)
    d = y(na+1) - ym;
    
    Y0 = ones(N,1)*d;
    for p=1:N
        for i=1:nb
            if( p-i <= -1 )
                Y0(p) = Y0(p) + b(i)*u((nb+1)+p-i); % u(k+p-i) -k=nb+1-> u(nb+1+p-i)
            else
                Y0(p) = Y0(p) + b(i)*u((nb+1)-1);   % u(k-1) -k=nb+1-> u(nb+1-1)
            end 
        end
        for i=1:na
            if( p-i <= 0 )
                Y0(p) = Y0(p) - a(i)*y((na+1)+p-i); % y(k+p-i) -k=na+1-> y(na+1+p-i)        
            else
                Y0(p) = Y0(p) - a(i)*Y0(p-i);       % Y0(k+p-i) -k=0-> Y0(p-i)
            end
        end 
    end
    
    % Classical form -- not used
    %dU = K*(Yzad-Y0); 
    %du = dU(1);
    
    du = Knu*(Yzad-Y0);
end

function u_y_check(a,b,u,y)


end