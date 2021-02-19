function [outRow,outCol,outMtx] =  writeSparseMM(rowH,colH,inMatrix,stubname)
    
    if nargin < 4 || isempty(stubname)
        error('Output file name needs to be specified');
    end

    if exist(stubname,'dir')
        dirname = stubname;
        stubname = 'out';
    elseif strcmp(stubname(end),'/')
        dirname = stubname;
        stubname = 'out';
        % [s,~,~] = mkdir(dirname);            
    else
        [dirname,stubname,suffix] = fileparts(stubname);
        stubname = [stubname suffix];
    end
          
    [D,N] = size(inMatrix);
    
    if D ~= length(rowH) 
        error('Row labels do no match number of rows in data matrix');
    end

    if N ~= length(colH) 
        error('Column labels do no match number of columns in data matrix');
    end
    
    if ~exist(dirname,'dir')                
        [s,~,~] = mkdir(dirname);            
    end      

    % Write col h
    outRow = [ dirname filesep stubname '.rowh.tsv' ];
    
    fh = fopen(outRow,'w');
    cellfun(@(x)fprintf(fh,'%s\n',x),rowH);  
    fclose(fh);
    
    outCol = [ dirname filesep stubname '.colh.tsv' ];
    fh = fopen(outCol,'w');
    cellfun(@(x)fprintf(fh,'%s\n',x),colH);  
    fclose(fh);
    
    % matrix.mtx
    [i,j,v] = find(inMatrix);
    zOut = [ i j v ]';    
    
    
    outMtx = [ dirname filesep stubname 'matrix.mtx' ];    
    fh = fopen(outMtx,'w');    
    fprintf(fh,'%%%%MatrixMarket matrix coordinate real general\n%%\n%d %d %d\n',D,N,length(i));    
    % arrayfun(@(x,y,z)fprintf(fh,'%d %d %d\n',x,y,z),i,j,v);
    fprintf(fh,'%d %d %d\n',zOut(:));
    fclose(fh);
end