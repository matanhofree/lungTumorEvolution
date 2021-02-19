function write_sge_xgboost_subtype(dataFile,clustFile,varFile,clVar,batchVar,outDir,outRun,inOpts)

    defaultOpts.ncores = 4;
    defaultOpts.isCategory = 1;
    defaultOpts.ugerMem = 20;
    defaultOpts.ugerTime = '24:00:00';
    defaultOpts.testFiles = 1;
    defaultOpts.dataMat = 'normTPM';
    
    defaultOpts.littleR = '/ahg/regevdata/users/mhofree/home_unix_user/R/conda_rhel7/R361/littler/bin/r';
    defaultOpts.sourceR = '$HOME/projects/cancer_SC/code/rfSignatureDiscovery/xgboostPredClass.R';            
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts)

    if ~exist(dataFile,'file')
        error('Data file does not exist');
    end
    
    if ~exist(clustFile,'file')
        error('Clsut file does not exist');
    end
    
%     if ~exist(varFile,'file') && ~isempty(varFile)
%         error('Var file does not exist');
%     end
    
    clRaw = fastTxtRead(clustFile,',');
    clV = data2simpleStruct(clRaw(2:end,:),clRaw(1,:))

    if ~isfield(clV,clVar)
        error('Variable does not exist in file');        
    end
            
    if ~isfield(clV,batchVar)
        error('Batch variable does not exist in file');        
    end
    
    tabFilter(clV.(clVar));
    
    cN = (luniq(clV.(clVar))-1)*(luniq(clV.(batchVar))+1);
    fprintf('Writing run file for %d jobs\n',cN);   
    
    
    fout = fopen(outRun,'w');
    fprintf(fout,'#!/bin/bash\n');
    fprintf(fout,'#\n');
    fprintf(fout,'#$ -cwd\n');
    fprintf(fout,'#$ -V\n');
    fprintf(fout,'#$ -q broad\n');
    fprintf(fout,'#$ -P regevlab\n');
    fprintf(fout,'#$ -e logs/\n');
    fprintf(fout,'#$ -o logs/\n');
    fprintf(fout,'#$ -pe smp %d\n',opts.ncores);
    fprintf(fout,'#$ -binding linear:%d\n',opts.ncores);
    fprintf(fout,'#$ -N xgb_%s_%s\n',clVar,batchVar);
    fprintf(fout,'#$ -l h_vmem=%dg\n',opts.ugerMem);
    fprintf(fout,'#$ -l h_rt=%s\n',opts.ugerTime);
    fprintf(fout,'#$ -t 1-%d\n',cN);
    fprintf(fout,'#$ -R y\n');
    fprintf(fout,'\n');
    fprintf(fout,'\n');
    fprintf(fout,'source /broad/software/scripts/useuse\n');
    fprintf(fout,'use miniconda\n');
    fprintf(fout,'source activate R361\n');
    fprintf(fout,'\n');
    fprintf(fout,'\n');
    if ~isempty(varFile)
        fprintf(fout,'VARGENES=%s\n',varFile);
    end
    fprintf(fout,'DATAFILE=%s\n',dataFile);
    fprintf(fout,'CLUSTFILE=%s\n',clustFile);
    fprintf(fout,'\n');
    fprintf(fout,'\n');   
    fprintf(fout,'OUTDIR=%s\n',outDir);
    fprintf(fout,'mkdir -p ${OUTDIR}\n');
    fprintf(fout,'OUTSTUB=%s/clP_%s_%s\n',outDir,clVar,batchVar);    
    fprintf(fout,'\n');       
    fprintf(fout,'\n');
    fprintf(fout,'\n');
    fprintf(fout,'LITTLER=%s\n',opts.littleR);
    fprintf(fout,'\n');
    fprintf(fout,'hostname\n');
    fprintf(fout,'date\n');
    fprintf(fout,'echo ${VARGENES}\n');
    fprintf(fout,'echo ${DATAFILE}\n');
    fprintf(fout,'echo ${CLUSTFILE}\n');
    fprintf(fout,'echo ${OUTSTUB}\n');
    fprintf(fout,'\n');       
    fprintf(fout,'\n');
    % fprintf(fout,'%s --vanilla --rtemp -e "source(''%s'');run_xgboost_classify(''${INDATA}'',''${OUTFILE}'',''${OPTSJSON}'')"\n',opts.littleR,opts.sourceR);
    fprintf(fout,'${LITTLER} --vanilla --rtemp -p <<_END_RSCRIPT_\n');
    fprintf(fout,'\n');       
    fprintf(fout,'source("%s")\n',opts.sourceR);
    fprintf(fout,'.libPaths()\n');
    fprintf(fout,'\n');
    fprintf(fout,'\n');
    fprintf(fout,'opts=list(dataMat="%s")\n',opts.dataMat);
    % xgboostPredClass("${DATAFILE}","${CLUSTFILE}",clustTableDim = "cNMFnetST12",outStub="${OUTSTUB}",runIdx=${SGE_TASK_ID},batchCV = "Kfold5",optsJson=opts)
    if ~isempty(varFile)    
        fprintf(fout,'xgboostPredClass("${DATAFILE}","${CLUSTFILE}",clustTableDim = "%s", geneFeatureFile = "${VARGENES}",outStub="${OUTSTUB}",runIdx=${SGE_TASK_ID},batchCV = "%s",optsJson=opts)\n',clVar,batchVar);
    else
        fprintf(fout,'xgboostPredClass("${DATAFILE}","${CLUSTFILE}",clustTableDim = "%s",outStub="${OUTSTUB}",runIdx=${SGE_TASK_ID},batchCV = "%s",optsJson=opts)\n',clVar,batchVar);
    end
    fprintf(fout,'\n');
    fprintf(fout,'\n');       
    fprintf(fout,'\n');
    fprintf(fout,'sessionInfo()\n');
    fprintf(fout,'_END_RSCRIPT_\n');
    fprintf(fout,'date\n');
    fprintf(fout,'\n');
    fclose(fout);   

end