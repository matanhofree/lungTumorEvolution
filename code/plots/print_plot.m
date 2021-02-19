function print_plot(cfig,fname,dirList,suffix,clobber,subDir,extOpts)

    if ~exist('suffix','var') || isempty(suffix)
        suffix = { '-dpng' };
    elseif ischar(suffix)
        suffix = { suffix };
    end
    
    if ~exist('dirList','var') || isempty(dirList)
        error('Specify output dir');
    elseif ischar(dirList)
        dirList = { dirList };
    end
    
    if ~exist('cfig','var') || isempty(cfig)
        cfig = { gcf };
    elseif ~iscell(cfig)
        cfig = {cfig};
    end
    
    if ~exist('clobber','var') || isempty(clobber)
        clobber = 0;
    end
    
    
    if ~exist('subDir','var') || isempty(subDir)
        if length(suffix) > 1
            subDir = 1;
        else
            subDir = 0;
        end
    end
    
    
    cN = length(cfig);
    if ~exist('fname','var') || isempty(fname)
        fname = cell(cN,1);
        % error('Need to name output');
    elseif ischar(fname)
        
%         fdir = regexprep(fname,'(.*)/[^/]+','$1');
        if cN > 1
%             if ~isempty(fdir)
%                 zsel = ~isemptycell(fdir);
%                 cellfun(@(x)mkdir(x),fdir,'unif',0);
%             end
            fname = mergeStringPair(fname,1:cN);
        else
            fname = { fname };
        end
    end
    
    if any(strgrep(fname,'/'))
        fdir = regexprep(fname,'(.+)/[^/]+','$1'); 
        fname =  regexprep(fname,'.+/([^/]+)','$1'); 
        
        dirList = mergeStringPair('%s/%s',dirList,fdir{1});
    else
        fdir = cell(cN,1);
    end
    
    
    for i = 1:length(dirList)        
        cDir = dirList{i};
        mkdir([ cDir ]);    
    end
    if subDir 
        for j = 2:length(suffix)
            for i = 1:length(dirList)        
                cDir = dirList{i};
                mkdir([ cDir '/' regexprep(suffix{j},'^-d','') ]);    
            end
        end
    end
    
    
    
    i = 1;
    for j = 1:length(suffix)
        cSuf = suffix{j};
        
        cDir = dirList{i};        
        if subDir && j > 1
            cDir = [ cDir '/' regexprep(suffix{j},'^-d','') ];
        end
        
        for zf = 1:length(cfig)

            cf = cfig{zf};
            cn = fname{zf};

            if isempty(cn)
                timeStamp = datestr(now,'yyyymmdd_HHMMSS');

                cName = matlab.lang.makeValidName(class(zf));
                outName = sprintf('%s/%s%s',cDir,cName,timeStamp); 
            else
                % timeStamp = datestr(now,'yyyymmdd_HHMMSS');
                outName = sprintf('%s/%s%s',cDir,cn); 
            end
            
            outFile = sprintf('%s.%s',outName,regexprep(cSuf,'-d',''));
            iName{j,zf} = regexprep(outFile,'\.epsc$','.eps');

            if clobber == 0 
                if exist(outFile,'file')
                    warning('%s file exists. Skipping',outFile);
                    continue;
                end
            end
            if strgrep(cSuf,'svg')
               cf.Renderer = 'painters'; 
            end
            fprintf('Saving plot: %s.%s\n',outFile,cSuf);
            print(cf,outName,cSuf);
            
            outFile = regexprep(outFile,'\.epsc$','.eps');

            % iName{j,zf} = outFile;
        end
    end
    for i = 2:length(dirList)       
        for j = 1:length(suffix)
            cDir = dirList{i};
            if subDir && j > 1
                cDir = [ cDir '/' regexprep(suffix{j},'^-d','') ];
            end
            cSuf = suffix{j};
            for zf = 1:length(cfig)

                cf = cfig{zf};
                cn = fname{zf};

                if isempty(cn)
                    timeStamp = datestr(now,'yyyymmdd_HHMMSS');

                    cName = matlab.lang.makeValidName(class(zf));
                    outName = sprintf('%s/%s%s',cDir,cName,timeStamp); 
                else
                    % timeStamp = datestr(now,'yyyymmdd_HHMMSS');
                    outName = sprintf('%s/%s%s',cDir,cn); 
                end

                outFile = sprintf('%s.%s',outName,regexprep(cSuf,'-d',''));
                outFile = regexprep(outFile,'\.epsc$','.eps');

                if clobber == 0 
                    if exist(outFile,'file')
                        warning('%s file exists. Skipping',outFile);
                        continue;
                    end
                end
                if strgrep(cSuf,'svg')
                   cf.Renderer = 'painters'; 
                end
                
                
                fprintf('Copying plot: %s\n',outFile);                                
                [statW,msg] = copyfile(iName{j,zf},outFile);
                if statW == 0 
                    fprintf('Failed to copy %s - %s\n',outFile,msg);
                end
                % print(cf,outName,cSuf);                    
            end
        end
    end
    
end