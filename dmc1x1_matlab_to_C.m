clearvars;
load('dmc1x1_001');

fprintf('const int N = %d;\n',N);
fprintf('const int Nu = %d;\n',Nu);
fprintf('const int D = %d;\n',D);
fprintf('const float Ke = %ff;\n',Ke);

% Mp
fprintf('const float Mp[N*(D-1)] = {\n');
for i=1:N
    for j=1:(D-1)
        fprintf('%+.6e,',Mp(i,j));
    end
    fprintf('\n');
end
fprintf('};\n');

% Ku
fprintf('const float Ku[D-1] = {');
for i=1:(D-1)
    fprintf('%+.6e,',Ku(i));    
end
fprintf('};\n');
