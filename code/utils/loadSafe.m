function outD = loadSafe(loadFile,varToKeep,varToIgnore)

    outD = [];
    if isstruct(loadFile)
        outD = loadFile;
        varToLoad = fieldnames(outD);        
        varToLoadRef = varToLoad;
        
        if nargin > 1 && ~isempty(varToKeep)        
            varToLoad = intersect(varToLoad,varToKeep);
        end

        if nargin > 2 && ~isempty(varToIgnore)
            varToLoad = setdiff(varToLoad,varToIgnore);            
        end

        if ~isequal(varToLoad,varToLoadRef)
            outD = structSelectField(outD,varToLoad);            
        end
    elseif exist(loadFile,'file')
        fileVariables = whos('-file',loadFile);

        varToLoad = { fileVariables.name };

        if nargin > 1 && ~isempty(varToKeep)        
            varToLoad = intersect(varToLoad,varToKeep);
        end

        if nargin > 2 && ~isempty(varToIgnore)
            varToLoad = setdiff(varToLoad,varToIgnore);            
        end

        fprintf('Loading the following variables from mat file: %s\n',loadFile);
        disp(varToLoad);
        outD = load(loadFile,varToLoad{:});
    else
        fprintf('No file or variable to load\n');
    end
end