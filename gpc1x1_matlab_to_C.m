clearvars;
load('gpc1x1_001');

fprintf('const int N = %d;\n',N);
fprintf('const int Nu = %d;\n',Nu);

na = length(a);
nb = length(b);
fprintf('const int na = %d;\n',na);
fprintf('const int nb = %d;\n',nb);

% a
fprintf('const float a[na] = {');
for i=1:na
    fprintf('%+.6e,',a(i));    
end
fprintf('};\n');
% b
fprintf('const float b[nb] = {');
for i=1:nb
    fprintf('%+.6e,',b(i));    
end
fprintf('};\n');

% Knu
fprintf('const float Knu[N] = {');
for i=1:N
    fprintf('%+.6e,',Knu(i));    
end
fprintf('};\n');
