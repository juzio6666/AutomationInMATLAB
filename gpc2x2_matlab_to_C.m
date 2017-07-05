% TODO !!!
stuff = '';
%clearvars;
%load('gpc2x2_001');

stuff = [stuff, sprintf('#define N %d\n',N)];
stuff = [stuff, sprintf('#define Nu %d\n',Nu)];

na = length(a);
nb = length(b);
stuff = [stuff, sprintf('#define na %d\n',na)];
stuff = [stuff, sprintf('#define nb %d\n',nb)];

stuff = [stuff, sprintf('#define nu %d\n',nu)];
stuff = [stuff, sprintf('#define ny %d\n',ny)];

stuff = [stuff, sprintf('const float b[ny][nu][nb] = {')];
for m=1:ny
    stuff = [stuff, '{'];
    for n=1:nu
        stuff = [stuff, '{'];
        for i=1:nb
            stuff = [stuff, sprintf('%+.6e',b(m,n,i))];
            if(i~=nb); stuff = [stuff, ',']; end
        end
        stuff = [stuff, '}'];
        if(n~=nu); stuff = [stuff, ',']; end
    end
    stuff = [stuff, '}'];
    if(m~=ny); stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];

stuff = [stuff, sprintf('const float a[ny][na] = {')];
for m=1:ny
    stuff = [stuff, '{'];
    for i=1:na
        stuff = [stuff, sprintf('%+.6e',a(m,i))];   
        if(i~=na); stuff = [stuff, ',']; end
    end     
    stuff = [stuff, '}'];
    if(m~=ny); stuff = [stuff, ',']; end
end 
stuff = [stuff, sprintf('};\n')];

% Knu
stuff = [stuff, sprintf('float Knu[nu*N*ny] = {')];
for n=1:nu
    for i=1:N*ny
        stuff = [stuff, sprintf('%+.6e',Knu(n,i))];  
        if(i~=N*ny || n~=nu); stuff = [stuff, ',']; end  
    end 
    if(n~=nu); stuff = [stuff, newline]; end
end
stuff = [stuff, sprintf('};\n')];
fprintf(stuff);
clipboard('copy', stuff);