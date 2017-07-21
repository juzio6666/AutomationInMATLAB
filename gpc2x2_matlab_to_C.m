stuff = '';

na = size(a,2);
nb = size(b,3);
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
            stuff = [stuff, sprintf('%+.6ef',b(m,n,i))];
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
        stuff = [stuff, sprintf('%+.6ef',a(m,i))];   
        if(i~=na); stuff = [stuff, ',']; end
    end     
    stuff = [stuff, '}'];
    if(m~=ny); stuff = [stuff, ',']; end
end 
stuff = [stuff, sprintf('};\n')];

% Kyzad
stuff = [stuff, sprintf('float Kyzad[nu][ny] = {')];
for n=1:nu
    stuff = [stuff, '{'];
    for i=1:ny
        stuff = [stuff, sprintf('%+.6ef',Kyzad(n,i))];  
        if(i~=ny); stuff = [stuff, ',']; end  
    end 
    stuff = [stuff, '}'];
    if(n~=nu); stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];

% Ky
stuff = [stuff, sprintf('float Ky[nu][ny][na+1] = {')];
for n=1:nu
    stuff = [stuff, '{'];
    for i=1:ny
        stuff = [stuff, '{'];
        for j=1:(na+1)
            stuff = [stuff, sprintf('%+.6ef',Ky(n,i,j))];  
            if(j~=(na+1)); stuff = [stuff, ',']; end
        end
        stuff = [stuff, '}'];
        if(i~=ny); stuff = [stuff, ',']; end
    end 
    stuff = [stuff, '}'];
    if(n~=nu); stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];
% Ky
stuff = [stuff, sprintf('float Ku[nu][nu][nb] = {')];
for n=1:nu
    stuff = [stuff, '{'];
    for i=1:nu
        stuff = [stuff, '{'];
        for j=1:nb
            stuff = [stuff, sprintf('%+.6ef',Ku(n,i,j))];  
            if(j~=nb); stuff = [stuff, ',']; end
        end
        stuff = [stuff, '}'];
        if(i~=nu); stuff = [stuff, ',']; end
    end 
    stuff = [stuff, '}'];
    if(n~=nu); stuff = [stuff, ',']; end
end
stuff = [stuff, sprintf('};\n')];
fprintf(stuff);
clipboard('copy', stuff);