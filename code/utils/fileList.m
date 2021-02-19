function [fList,fPath,fName] = fileList(dirName,suffix)
    
% zbuf = ls(dirName);
% flist = textscan(zbuf,'%s');
% flist = flist{:};
% 
% for i = 1:length(flist)
%     [fPath{i},zn,zext]= fileparts(flist{i});
%     fName{i} = [ zn zext ];
% end
fList = [];
fPath = [];
fName = [];
[inPath] = fileparts(dirName);

if ~exist('suffix','var')
    suffix = '*';
    zbuf = dir([dirName suffix]);
elseif isempty(suffix)   
    zbuf = dir(dirName);
else
    zbuf = dir([dirName suffix]);
end


fName = { zbuf(~[zbuf.isdir]).name }';
fPath = { zbuf(~[zbuf.isdir]).folder }';

fList = cellfun(@(x,y)sprintf('%s%s%s',x,filesep,y),fPath,fName,'uniformoutput',0);


