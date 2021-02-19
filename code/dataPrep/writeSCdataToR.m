function writeSCdataToR(countData,outDir,inName,skipVars,geneID,sampleID,inOpts)
    
    defaultOpts.defName = 'countData';
    defaultOpts.gzip = 1;
    defaultOpts.sampleKey = 'sampleID';
    defaultOpts.geneKey = 'ensgID';
    defaultOpts.makeGeneKeyUniq = 1;
    defaultOpts.structSep = '_dZ_';
    defaultOpts.sparseSep = '_dSp_';
    defaultOpts.denseSep = '_dDD_';

    defaultOpts.saveExtraAnnot = 1;
    
    defaultOpts.writeSparse = 'HDF5'; % 'MM' or 'HDF5'
    
    defaultOpts.maxDenseCol = 20;
    defaultOpts.removeZeros = 1;
    defaultOpts.inSub = 0;
    defaultOpts.writeDenseMat = 0;
    
    defaultOpts.skipData = 0;
    
    defaultOpts.tableAsTxt = 1;
    defaultOpts.denseAddSmpHeader = 1;
    defaultOpts.delim = ',';
    defaultOpts.suffix = '.csv';
                
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);      
    
    if ~exist('skipVars','var') || isempty(skipVars)
        skipVars = {};        
    end
    
    if ~exist('inName','var') || isempty(inName)
        inName = inputname(1);     
    end
    
    if isempty(inName)       
        inName = opts.defName;        
    end           
    
    if ~exist('sampleID','var') || isempty(sampleID)        
        sampleID = countData.(opts.sampleKey);
    end
    if ~exist('geneID','var') || isempty(geneID)        
        geneID = countData.(opts.geneKey);
    end
    
    N = length(sampleID);    
    D = length(geneID);
    
    if opts.makeGeneKeyUniq && D > luniq(geneID)
        geneID = matlab.lang.makeUniqueStrings(geneID);
    end
    
    skipVars{end+1} = opts.sampleKey;
    % skipVars{end+1} = opts.geneKey;
        
    if ~exist(outDir,'dir') 
        [~,~,~] = mkdir(outDir);
    end
    
    if isempty(opts.suffix)
        if strcmp(opts.delim,'\t')
            opts.suffix = '.tsv';
        else 
            opts.suffix = '.csv';
        end
    end
    
    fnames = fieldnames(countData); 
    
    sampleRowHeader = cell2table(sampleID,'VariableNames',{ opts.sampleKey });
    geneRowHeader = cell2table(geneID,'VariableNames',{ opts.geneKey });
    
    sampleColHeader = matlab.lang.makeValidName(sampleID);
    geneColHeader = matlab.lang.makeValidName(geneID);
    
    
    smpTable.sampleID = sampleID;
    geneTable.geneID = geneID;    
    cSmp = 1;
    cGene = 1;
     
    
    for i = 1:length(fnames)
        cname = fnames{i};   
        
        if ismember(cname,skipVars)            
            fprintf('Skipping: %s (inList)\n',cname);
            continue;
        end
                
        cMat = countData.(cname);
        [zD,zN] = size(cMat);                        
        
        if isstruct(cMat) 
            opts.inSub = 1;
            subName = [ outDir '/' inName opts.structSep cname ];
            writeSCdataToR(cMat,subName,cname,skipVars,geneID,sampleID,opts);
        elseif issparse(cMat) && (zD == D) && (zN == N)
            
            subName = [ outDir '/' inName opts.sparseSep cname ];
            
            if opts.skipData 
                fprintf('Skipping writing of data matrix (%s)\n',cname);                    
                continue;
            end
            
            switch opts.writeSparse                
                case 'MM'
                    fprintf('Writing %s in MM format\n',cname);                    
                    [outRow,outCol,outMtx] = writeSparseMM(geneID,sampleID,cMat,subName);             

                    if opts.removeZeros
                        if ismac()
                            system([ '/usr/local/bin/gsed -i ''s/0.000000/0/g'' ' outMtx ]);
                        else
                            system([ '/bin/sed -i ''s/0.000000/0/g'' ' outMtx ]);
                        end
                    end
                    if opts.gzip
                        gzip(outMtx);
                    end
                case 'HDF5'
                    fprintf('Writing %s in hdf5 (matlab) format\n',cname);                    
                    outS = [];
                    [ outS.i, outS.j, outS.v ] = find(cMat);
                    save(subName,'-v7.3','-struct','outS');                    
                    
%                 case 'MAT'
%                     fprintf('Writing %s in Matlab format\n',cname);                    
%                     outS = [];
%                     [ outS.i, outS.j, outS.v ] = find(cMat);
%                     
%                     outS.geneID = geneID;
%                     outS.sampleID = sampleID;
%                     
%                     save(subName,'-v7.3','-struct','outS');
                    
                otherwise
                    fprintf('Skipping: %s. Choose a sparse format to use\n',cname);
                    
            end
    
        else
            % Not sparse             
            if min([zD zN]) == 1 % Vector
                if zD == D
                    geneTable.(cname) = cMat;
                    cGene = cGene + 1;
                elseif zD == N
                    smpTable.(cname) = cMat;
                    cSmp = cSmp + 1;
                elseif zN == N 
                    smpTable.(cname) = cMat';
                    cSmp = cSmp + 1;
                elseif zN == D
                    geneTable.(cname) = cMat';
                    cGene = cGene + 1;
                end
            else % Matrix
                if opts.saveExtraAnnot
                    
                    % subName = [ outDir '/' inName '_' cname '.csv' ];

                    if ~isempty(opts.maxDenseCol) && min([zD zN]) < opts.maxDenseCol
                        
                        subName = [ outDir '/' inName opts.denseSep cname opts.suffix ];

                        fprintf('Writing table %s -- %s\n',cname,subName);
                        writetable(array2table(cMat),subName);       
                    elseif opts.writeDenseMat
                        if ~istable(cMat)
                            if iscell(cMat)
                                cMat = cell2table(cMat);
                            else
                                cMat = array2table(cMat);
                            end  
                            if opts.denseAddSmpHeader
                                if zN == N 
                                    cMat.Properties.VariableNames = sampleColHeader;
                                elseif zN == D
                                    cMat.Properties.VariableNames = geneColHeader;
                                end 
                            end
                        end
                        
                        if opts.denseAddSmpHeader
                            if zD == N 
                                cMat = [ sampleRowHeader cMat ];
                            elseif zD == D
                                cMat = [ geneRowHeader cMat ];
                            end 
                        end
                                                  
                        if opts.tableAsTxt == 1

                            subName = [ outDir '/' inName opts.denseSep cname opts.suffix ];
                            
                            fprintf('Writing table (%dx%d) as text %s -- %s\n',zD,zN,cname,subName);
                            writetable(cMat,subName,'Delimiter',opts.delim);                            
                        else 
                            subName = matlab.lang.makeValidName( [ inName cname ] );
                            outS.(subName) = cMat;
                            fprintf('Writing table (%dx%d) as h5 %s -- %s\n',zD,zN,cname,subName);

                            subName = [ outDir '/' inName opts.denseSep cname '.h5' ];

                            fprintf('Writing h5 table %s -- %s\n',cname,subName);
                            save(subName,'-v7.3','-struct','outS');  
                        end
                        
                    else
                        fprintf('Skipping variable: %s (%d x %d). Number of columns is greater than maxDenseCol.\n',cname,zD,zN)                                
                    end
                else
                    fprintf('Skipping variable: %s\n..\n\tSet saveExtraAnnot to true to save.\n',cname);
                end                                                                                                            
            end
        end
    end 
    
    if opts.inSub == 0 || cSmp > 1 
        subName = [ outDir '/' inName '_smpTable' opts.suffix ];     
        writetable(struct2table(smpTable),subName,'Delimiter',opts.delim);  
    end
        
    if opts.inSub == 0 || cGene > 1
        subName = [ outDir '/' inName '_geneTable' opts.suffix ];
        writetable(struct2table(geneTable),subName,'Delimiter',opts.delim);  
    end
end


