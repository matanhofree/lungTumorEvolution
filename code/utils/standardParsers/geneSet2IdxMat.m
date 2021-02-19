function [idxGeneKey,idxMat] = geneSet2IdxMat(setGeneList)    

    allGeneList = concatCell_to_cell(setGeneList);
    idxGeneKey = unique(allGeneList); 
    lenN = length(setGeneList);
    idxMat = false(lenN,length(idxGeneKey));    
    % keyPosHash = containers.Map(idxGeneKey ,1:length(idxGeneKey));    
   
    %%
    tic;
    zposAll = fasthash(idxGeneKey,allGeneList);      
    toc; 
    %%
    gsetLenAll = cellfun(@(x)length(x),setGeneList);
    gsetLenVect = arrayfun(@(x,y)repmat(x,y,1),(1:lenN)',gsetLenAll,'uniformoutput',0);
    gsetLenVect = vertcat(gsetLenVect{:});
    
    idxMat = sparse(gsetLenVect,zposAll,1,lenN,length(idxGeneKey));
      
%     progressbar();
%     zcnt = 1;
%     for i = 1:lenN        
%         % zpos = fasthash(idxGeneKey,setGeneList{i});        
%         zpos = zposAll(zcnt+)
%         idxMat(i,zpos) = 1;
%         progressbar(i/lenN);
%     end