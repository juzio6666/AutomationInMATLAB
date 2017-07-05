% TODO !!!
stuff = '';
%clearvars;
%load('gpc2x2_001');

stuff = [stuff, sprintf('#define N %d\n',N)];
stuff = [stuff, sprintf('#define Nu %d\n',Nu)];

stuff = [stuff, sprintf('#define nu %d\n',nu)];
stuff = [stuff, sprintf('#define ny %d\n',ny)];

% Ku
stuff = [stuff, sprintf('float Ku[nu*N*ny] = {')];
for n=1:nu
    for i=1:N*ny
        stuff = [stuff, sprintf('%+.6e',Knu(n,i))];  
        if(i~=N*ny || n~=nu); stuff = [stuff, ',']; end  
    end 
    if(n~=nu); stuff = [stuff, newline]; end
end
stuff = [stuff, sprintf('};\n')];

clipboard('copy', stuff);