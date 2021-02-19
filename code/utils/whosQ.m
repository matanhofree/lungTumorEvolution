function sTable=whosQ(varargin)
%Whos (Upper case W) is  a customized version of whos (lower case w).
%
%The differences are 
%
% (1) "Whos" will compute/display memory in kilobytes whereas "whos" displays in bytes.
%
% (2) "Whos" will always display all array dimensions whereas "whos" will
%     not display array dimensions for 4th and higher dimensional arrays.
%     
%
%EXAMPLE: Given arrays
%
%  A=rand(10,20,30,40); 
%  B=rand(5,10,15,20,25,'single')*i;
%
%using "whos" displays the following:
%
%>> whos
%   Name      Size              Bytes  Class     Attributes
% 
%   A         4-D             1920000  double              
%   B         5-D             3000000  single    complex   
%
%
%whereas "Whos" will display the following
%
%>> Whos
%   Name   Size              Kilobytes     Class    Attributes
%                                                             
%   A      10x20x30x40            1875     double             
%   B      5x10x15x20x25          2930     single   complex   
%
%
%An unfortunate limitation is that Whos relies on EVALIN and so will not work
%correctly in a workspace context reached using DBUP and DBDOWN.
%I elaborate on this somewhat in this Newsgroup thread:
%
%      http://www.mathworks.com/matlabcentral/newsreader/view_thread/303357
%
%Hopefully, TMW will provide a way around this eventually.
%
%by Matt Jacobson
%  
%Copyright, Xoran Technologies, Inc. 2011
  
if isempty(varargin)
    cmdstr='whos;';
    s=evalin('caller', cmdstr) ;
    %%

    % SS = cellfun(@(x)(  s.(x) ' ),fieldnames(s),'uniformoutput',0)'
    % 
    SS = struct2cell(s)';
else
    for i = 1:length(varargin)
        if isstruct(varargin{i})
            zvar = inputname(i);
            if isempty(zvar)
               zvar = sprintf('unlabeled_%d',i); 
            end            
            eval(sprintf('%s = varargin{%d}',zvar,i));
            eval(sprintf('structunpack(%s)',zvar));                        
        else 
            error('Not implemented');
        end
    end
    
    s = whos();
    SS = struct2cell(s)';
end

%%
znames=fieldnames(s);
znames(5:end) = regexprep(znames(5:end),'(.*)','is_$1');

% drop nesting
znames(8) = [];
SS(:,8) = [];
%%
% Fix dim 
dimV = SS(:,2);
isSingle = cellfun(@(x)all(x==1),dimV);
isTensor = cellfun(@(x)length(x)>2,dimV);
is2D = cellfun(@(x)size(x,2)==2,dimV);
%
dimStr = cell(size(dimV));
dimStr(:) = { '?' };
dimStr(isSingle) = {'1'};
dimStr(is2D) = cellfun(@(x)sprintf('%dx%d',x),dimV(is2D),'uniformoutput',0);
dimStr(isTensor) = cellfun(@(x)[ sprintf('%d',x(1)) sprintf('x%d',x(2:end)) ],dimV(isTensor),'uniformoutput',0);

%%
SS(:,2) = dimStr;

%% Fix size string

varSize = cell2mat(SS(:,3));
%%
sizeStr = cell(size(varSize));
for ii = 1:size(varSize,1)
    cur_size = varSize(ii);
    if cur_size < 1025
        sizeStr{ii} = sprintf('%dB', cur_size);
    elseif cur_size < 1024^2+1
        sizeStr{ii} = sprintf('%dKB', round(cur_size/1024));
    elseif cur_size < 1024^3+1
        sizeStr{ii} = sprintf('%dMB', round(cur_size/1024^2));
    else
        sizeStr{ii} = sprintf('%dGB', round(cur_size/1024^3));
    end
end
%%
SS(:,3) = sizeStr;

[~,zidx] = sort(varSize);

SS = SS(zidx,:);

%%
sTable = cell2table(SS,'variablenames',znames);
% if nargout < 1
%     disp(sTable);
% end
