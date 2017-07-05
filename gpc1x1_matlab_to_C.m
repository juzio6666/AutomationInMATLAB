stuff = '';

stuff = [stuff, sprintf('const int N = %d;\n',N)];
stuff = [stuff, sprintf('//const int Nu = %d;\n',Nu)];

stuff = [stuff, sprintf('const int na = %d;\n',na)];
stuff = [stuff, sprintf('const int nb = %d;\n',nb)];

% a
stuff = [stuff, sprintf('const float a[na] = {')];
for i=1:na
    stuff = [stuff, sprintf('%+.6e,',a(i))];    
end
stuff = [stuff, sprintf('};\n')];
% b
stuff = [stuff, sprintf('const float b[nb] = {')];
for i=1:nb
    stuff = [stuff, sprintf('%+.6e,',b(i))];    
end
stuff = [stuff, sprintf('};\n')];

% Knu
stuff = [stuff, sprintf('const float Knu[N] = {')];
for i=1:N
    stuff = [stuff, sprintf('%+.6e,',Knu(i))];    
end
stuff = [stuff, sprintf('};\n')];

fprintf(stuff);
clipboard('copy', stuff);