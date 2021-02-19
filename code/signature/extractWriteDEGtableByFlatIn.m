function [outDeg,summaryTable] = extractWriteDEGtableByFlatIn(degTable,geneID,geneFilter,geneOrder,outPath,inOpts)

    
    defaultOpts.doWrite = 1;
    defaultOpts.doSummary = 1;
    defaultOpts.summaryGenes = 100;
    defaultOpts.filterSubset = 100;   
    defaultOpts.singleXLS = 1;   
    
    defaultOpts.clustNames = 'outClustNames';    
    defaultOpts.clNamesMap = [];
    defaultOpts.fieldList = { 'pval' 'fdr' 'meanExp' 'logR' 'expT' 'rsZ' 'auc' 'subMinAuc' 'subMinIdx' 'maxLogR' };        
    defaultOpts.fieldListRename = []; 
    defaultOpts.resolveIdx = [];
    defaultOpts.dropEmptyP = 1;
     
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end    
    clear defaultOpts;
  
    
    if ~exist('geneFilter','var')    
        geneFilter = [];
    end

    if ~exist('geneOrder','var')    
        geneOrder = [];         
    end

    if ~exist('outPath','var')    
        outPath = []; 
        opts.doWrite = 0;
        fprintf('Note - not writing output files outPath not specfied');
    end

    if opts.doWrite
        if opts.singleXLS 
            if ~isempty(opts.filterSubset) & ~isinf(opts.filterSubset) & ~isempty(geneOrder)
                outFile = sprintf('%s_top%d.xlsx',outPath,opts.filterSubset);
            else
                outFile = sprintf('%s_all.xlsx',outPath);
            end
        else            
            if ~isempty(opts.filterSubset) & ~isinf(opts.filterSubset) & ~isempty(geneOrder)
                outFile = sprintf('%s_top%d',outPath,opts.filterSubset);
            else
                outFile = sprintf('%s_all',outPath);
            end            
        end        
    end
    
    if ~isempty(opts.fieldList)
        fieldList = opts.fieldList;
        fieldList(~ismember(fieldList,fieldnames(degTable))) = [];
        degTableSub = structSelectField(degTable,fieldList);
   
        degTableSub = orderfields(degTableSub,fieldList);        
    else                
        degTableSub = degTable;   
        fieldList = fieldnames(degTableSub);
    end    
        
    if ~isempty(geneID)
        fieldList = [ 'geneID'; fieldList(:) ];
        nG = length(geneID);
    end
    
    
    if ~isempty(geneFilter)
        fieldList = [ fieldList(:); 'geneFilter' ];
        degTableSub.geneFilter = geneFilter;
    end
    
    if isnumeric(degTable.(opts.clustNames))
        clustNames = matlab.lang.makeValidName(num2cellstr(degTable.(opts.clustNames)));        
    else
        clustNames = matlab.lang.makeValidName(degTable.(opts.clustNames));        
    end    
    
    nC = length(clustNames);        
    
    
    if nC == 2
        zfnames = fieldnames(degTableSub);
        emptyField = structfun(@(x)isempty(x),degTableSub);
        degTableSub = rmfield(degTableSub,zfnames(emptyField));
        
        fieldList(ismember(fieldList,zfnames(emptyField))) = [];
    end
    
    if opts.doSummary
        summaryTable = cell(opts.summaryGenes,nC);
    else
        summaryTable = [];
    end
    
    
    for i = 1:nC        
        %%
        subDeg = structfun(@(x)x(:,i),degTableSub,'uniformoutput',0);
               
        if ~isempty(geneID)
            subDeg.geneID = geneID;
        end

        subDeg = orderfields(subDeg,fieldList);
        
        geneO = [];
        if ~isempty(geneFilter) && ~isempty(opts.filterSubset)
            zSelG = geneFilter(:,i)>0;
            
            fprintf('C%d) found %d genes\n',i,sum(zSelG));
            if sum(zSelG) == 0
                fprintf('Skipping');
                continue;
            end
            
            subDeg = structSubSelectMat(subDeg,zSelG);            
            if ~isempty(geneOrder) 
                [~,geneO] = sort(geneOrder(zSelG,i),'descend','MissingPlacement','last');
            end            
        else
            if ~isempty(geneOrder) 
                [~,geneO] = sort(geneOrder(:,i),'descend','MissingPlacement','last');
            end                        
        end             
        
        if ~isempty(opts.resolveIdx) && isfield(subDeg,opts.resolveIdx)
            rIdx = subDeg.(opts.resolveIdx);
            cIdxRes = cell(length(rIdx),1); 
            cIdxRes(:) = { '' };
            
            cIdxRes(~isnan(rIdx)) = clustNames(rIdx(~isnan(rIdx)));
            subDeg.(opts.resolveIdx) = cIdxRes;
            fprintf('Sim to %s\n:',clustNames{i});
            tabFilter(cIdxRes);
        end               
        
        subDeg = struct2table(subDeg);
        
        if ~isempty(geneO)
            subDeg = subDeg(geneO,:);
        end        
        %%        
        subDeg = table2struct(subDeg,'toScalar',1);
        
        if opts.doSummary            
            sumGeneList = subDeg.geneID;
            if ~isempty(geneFilter)                   
                sumGeneList(subDeg.geneFilter ~= 1) = [];            
            end
            nLeft = size(sumGeneList,1);             
            summaryTable(1:min(opts.summaryGenes,nLeft),i) = sumGeneList(1:min(opts.summaryGenes,nLeft));
        end
        
        if opts.dropEmptyP
            zSelE = ~isnan(subDeg.pval);
            fprintf('Dropping empty P (%d)\n',sum(~zSelE));
            subDeg = structSubSelectMat(subDeg,zSelE);              
        end
        %%
        
        subDeg = struct2table(subDeg);
        
        if ~isempty(opts.fieldListRename)
           subDeg.Properties.VariableNames = opts.fieldListRename; 
        end
        
        nLeft = size(subDeg,1);                
        if ~isempty(opts.filterSubset)
             subDeg = subDeg(1:min(opts.filterSubset,nLeft),1:end-1);               
        end
        outDeg.(clustNames{i}) = subDeg;
        
        if ~isempty(outPath) && opts.doWrite
            if opts.singleXLS
                writetable(outDeg.(clustNames{i}),outFile,'Sheet',clustNames{i}(1:min(31,length(clustNames{i}))));             
            else
                outFileS = sprintf('%s_%s.csv',outFile,clustNames{i});
                writetable(outDeg.(clustNames{i}),outFileS);
            end                        
        end
        
    end
    
    if ~isempty(summaryTable)
        summaryTable = cell2table(summaryTable,'variableNames',clustNames);
    
    
        if ~isempty(outPath) && opts.doWrite
            if opts.singleXLS
                writetable(summaryTable,outFile,'Sheet',nC+1);             
            else
                outFileS = sprintf('%s_summary.csv',outFile);
                writetable(summaryTable,outFileS);
            end                        
        end
    end
end