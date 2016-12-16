function print_float_matrix(declaration, X)
    xname = sprintf('%s = {',declaration);
    fprintf(xname);
    for x=1:size(X,1)
        for y=1:size(X,2)
            if(y==1 && x==1); fprintf('');
            elseif(y==1); fprintf(',\n%s',(xname*0+1)*' ');
            else fprintf(',');
            end
            fprintf('%+2.10f',X(x,y));
        end
    end
    fprintf('};\n');
end