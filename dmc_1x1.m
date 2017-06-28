function du = dmc_1x1(S,N,Nu,Lambda,Psi,dUp,Y,Yzad)
%DMC_1X1 Dynamic Matrix Control.
%   du = DMC_1X1(S,N,Nu,Lambda,Psi,dUp,Y,Yzad) denotes new control value du 
%   for a process defined by step response S, by solving following problem:
%
%      min (Yzad-Y)'*Psi*(Yzad-Y) + du'*Lambda*du
%      du
%
%   Notice that there are no constraints imposed.
%
%   S is a vector of step response coefficients. Its size denotes the
%   dynamic horison D = length(S).
%   N is a scalar denoting prediction horison. 
%   Nu is a scalar denoting control horison.
%   Lambda and Psi are parameters used to change the influence of different
%   parts of minimised cost function. Lambda is a diagonal matrix with size
%   Nu x Nu.
%   dUp is a vector of past increments of control value. It is as follows:
%      
%      dUp = [du(k-1), du(k-2), ... , du(k-D+1)]'
%      
%   y is a scalar denoting a value of current measurement.
%   Yzad is a vector of length N. It denotes set point values over the
%   prediction horison. If Yzad is a scalar yzad, it is interpreted as:
%
%      Yzad = [yzad, yzad, ... , yzad]'
%
%   otherwise it should be in form of:
%
%      Yzad = [yzad(k+1), yzad(k+2), ... , yzad(k+N)]'


%% Here should be all of the correctness checking.
% S check
if(~isvector(S)); error('S should be a vector of any positive length');
end

% Defining the D
D = length(S);

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
% dUp check
if(isvector(dUp))
    if(length(dUp) ~= (D-1)); error('Vector dUp should have a size of D-1=%d instead of %d',D-1,length(dUp));
    elseif(size(dUp,1)~=(D-1)); dUp = dUp';
    end
else
    error('dUp is not a vector');
end

% Y check
if(isscalar(Y)); Y = ones(N,1)*Y;
elseif(isvector(Y))
    if(length(Y) ~= N); error('Wrong length of Y vector -- got %d expected N=%d',length(Y),N);
    else
        if(size(Y,1)==1); Y=Y';
        end
    end
else
    error('Y should be a vector of length N=%d',N);
end       

% Yzad check
if(isscalar(Yzad)); Yzad=ones(N,1)*Yzad;
elseif(isvector(Yzad))
    if(length(Yzad) ~= N); error('Wrong length of Yzad vector -- got %d expected N=%d',length(Yzad),N);
    elseif(size(Yzad,1) == 1); Yzad=Yzad'; 
    else % Yzad is OK as is
    end
else; error('Yzad should be a 1xN = 1x%d vector',N); 
end

%% Matrices, vectors and scalars preparation
dispall = 0;
if(dispall); N, Nu, S, Psi, Lambda, Yzad, Y, dUp, end

% Memmory assignment
M = zeros(N,Nu);
Mp = zeros(N,D-1);
Y0 = zeros(N,1);

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

%% Algorithm itself
Y0 = Y+Mp*dUp;
K = (M'*Psi*M+Lambda)^(-1)*M';
dU = K*(Yzad-Y0);
if(dispall); M,Mp,Y,Y0,K,dU, end
du = dU(1);
end
