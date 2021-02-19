function edgeList = ot_to_ancestorGraph_splitSummary(otMat,stageSet,clustV,stageList,stageV,inOpts)

    defaultOpts.normOut = 1;
    defaultOpts.normOutDir = 2;
    defaultOpts.minFrac = 0.1;
    defaultOpts.plotAll = 0;
    % defaultOpts.propTotal = 1;
    defaultOpts.connectMST = 2;
    defaultOpts.minClusterSize = 10;

    defaultOpts.postNormFull = 1;

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;


    sumRows = @(zdata,types,typeSelect)cell2mat(arrayfun(@(x)sum(zdata(types==x,:),1),(typeSelect),'uniformoutput',0));
    sumCols = @(zdata,types,typeSelect)cell2mat(arrayfun(@(x)sum(zdata(:,types==x),2),(typeSelect),'uniformoutput',0)');
    toCol = @(x)x(:);

    sourceList = [];
    targetList = [];
    valList = [];
    outSumPrev = [];

    for zi = 1:length(stageList)

        stageA = stageList{zi}(1);
        stageB = stageList{zi}(2);

        zAfull = ismember(stageV,stageA);
        zBfull = ismember(stageV,stageB);

        %%
        if iscell(stageA{1})
            stageA = stageA{1};
        end

        if iscell(stageB{1})
            stageB = stageB{1};
        end
        %%

        outTransport = otMat{zi}; % (zAfull,zBfull);

        zAclust = clustV(zAfull);
        zBclust = clustV(zBfull);
        %%

        if ~isempty(opts.minClusterSize)
            [namesA,~,~,cntA] = uniquec(zAclust);
            [namesB,~,~,cntB] = uniquec(zBclust);

            cselectA = sort(namesA(cntA>opts.minClusterSize));
            cselectB = sort(namesB(cntB>opts.minClusterSize));
            cselectA = cselectA(:);
            cselectB = cselectB(:);
        else
            cselectA = unique(zAclust);
            cselectB = unique(zBclust);
        end
        %%
        if isempty(cselectA) || isempty(cselectB)
            continue;
        end


        % aClustLabel = arrayfun(@(x)sprintf('%s_%02d',stageA{stA},x),cselectA,'uniformoutput',0);
        % bClustLabel = arrayfun(@(x)sprintf('%s_%02d',stageB{stB},x),cselectB,'uniformoutput',0);
        %%
        aClustLabel = arrayfun(@(x)sprintf('%s_cl%02d',stageA{1},x),cselectA,'uniformoutput',0);
        bClustLabel = arrayfun(@(x)sprintf('%s_cl%02d',stageB{1},x),cselectB,'uniformoutput',0);
        %%

        outSum = sumRows(sumCols(outTransport,zBclust,cselectB),zAclust,cselectA);
        %%
        if opts.normOut
            outSum = bsxfun(@times,outSum,sum(outSum,opts.normOutDir).^-1);
            outSumNorm = outSum;
        else
            outSumNorm = bsxfun(@times,outSum,sum(outSum,opts.normOutDir).^-1);
        end
        %%
        if ~isempty(opts.minFrac)

            [~,ziCol] = max(outSumNorm,[],2);
            % [~,ziRow] = max(outSumNorm,[],1);

            % connectMST = full(sparse(1:length(ziCol),ziCol,1,size(outSumNorm,1),size(outSumNorm,2))) | full(sparse(ziRow,1:length(ziRow),1,size(outSumNorm,1),size(outSumNorm,2)));
            dropEdges = outSumNorm<opts.minFrac;

            if opts.connectMST
                % connectMST = (full(sparse(ziRow,1:length(ziRow),1,size(outSumNorm,1),size(outSumNorm,2))));
                connectMST = (full(sparse(1:length(ziCol),ziCol,1,size(outSumNorm,1),size(outSumNorm,2))));
            end

            if opts.connectMST == 1
                dropEdges = dropEdges & connectMST == 0;
            end

            outSum(dropEdges) = 0;

            if opts.connectMST == 2
                outSum(connectMST & outSum == 0) = -1;
            end
        end

        %%
        if opts.plotAll
            outSum(outSum == 0) = -1;
        end
        if opts.postNormFull
            outSum = outSum/sum(outSum(:));
        end

        [i,j,v] = find(outSum);
        % if opts.plotAll
        v(v == -1) = 0.000001;
        % end


        sourceList = [sourceList; toCol(aClustLabel(i))];
        targetList = [targetList; toCol(bClustLabel(j))];


        valList = [valList; toCol(v)];

        % outSumPrev = outSum;
        %             end
        %         end
    end

    edgeList.sourceList = sourceList;
    edgeList.targetList = targetList;
    edgeList.valList = valList;

end