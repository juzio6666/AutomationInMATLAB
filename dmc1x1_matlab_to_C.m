stuff = '';

stuff = [stuff, sprintf('#define nu 1\n')];
stuff = [stuff, sprintf('#define ny 1\n')];
stuff = [stuff, sprintf('#define D %d\n',D)];

stuff = [stuff, sprintf('float Ke[nu*ny] = {%+.6ef};\n',Ke)];
% Ku
stuff = [stuff, sprintf('float Ku[nu*(D-1)*nu] = {')];
for i=1:(D-1)
    stuff = [stuff, sprintf('%+.6ef',Ku(i))];    
    if(i~=(D-1)) stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];

fprintf(stuff);
clipboard('copy', stuff);