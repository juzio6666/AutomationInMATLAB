stuff = '';

stuff = [stuff, sprintf('#define nu %d\n',nu)];
stuff = [stuff, sprintf('#define ny %d\n',ny)];

stuff = [stuff, sprintf('#define D %d\n',D)];

% Ku
stuff = [stuff, sprintf('float Ku[nu*(D-1)*nu] = {')];
for n=1:nu
    for i=1:(D-1)*nu
        stuff = [stuff, sprintf('%+.6e',Ku(n,i))];  
        if(i~=(D-1)*nu || n~=nu); stuff = [stuff, ',']; end  
    end 
    if(n~=nu); stuff = [stuff, newline]; end
end
stuff = [stuff, sprintf('};\n')];

% Ke
stuff = [stuff, sprintf('float Ke[nu*ny] = {')];
for n=1:nu
    for m=1:ny
        stuff = [stuff, sprintf('%+.6e',Ke(n,m))];  
        if(i~=ny || n~=nu); stuff = [stuff, ',']; end  
    end
    if(n~=nu); stuff = [stuff, newline]; end
end
stuff = [stuff, sprintf('};\n')];

fprintf(stuff);
clipboard('copy', stuff);