stuff = '';

na = length(a);
nb = length(b);
stuff = [stuff, sprintf('#define na %d\n',na)];
stuff = [stuff, sprintf('#define nb %d\n',nb)];

% a
stuff = [stuff, sprintf('const float a[na] = {')];
for i=1:na
    stuff = [stuff, sprintf('%+.6ef,',a(i))];    
end
stuff = [stuff, sprintf('};\n')];
% b
stuff = [stuff, sprintf('const float b[nb] = {')];
for i=1:nb
    stuff = [stuff, sprintf('%+.6ef,',b(i))];    
end
stuff = [stuff, sprintf('};\n')];

% Kyzad
stuff = [stuff, sprintf('const float Kyzad = %+.6ef;\n', Kyzad)];

% Ku
stuff = [stuff, sprintf('const float Ku[nb] = {')];
for i=1:nb
    stuff = [stuff, sprintf('%+.6ef',Ku(i))];    
    if(i~=nb) stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];


% Ku
stuff = [stuff, sprintf('const float Ky[na+1] = {')];
for i=1:(na+1)
    stuff = [stuff, sprintf('%+.6ef',Ky(i))];    
    if(i~=(na+1)) stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];

fprintf(stuff);
clipboard('copy', stuff);