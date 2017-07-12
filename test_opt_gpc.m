%function J = test_opt_gpc(Ku, Ky)
    global u y kp kk na nb Y0temp Y0new Y0temp2
%    Y0new = zeros(kk,1);
    for k=kp:kk
        X(:,k)=-[u(k-(1:nb));y(k-(0:na))];
%        Y0new(k) = - Ku*u(k-(1:nb)) - Ky*y(k-(0:na));
    end
    A=X'\Y0temp2;
    
%    J = Y0new-Y0temp2;
%    J = J'*J;
%end
