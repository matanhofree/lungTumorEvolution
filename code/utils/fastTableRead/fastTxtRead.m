function  [zdata,zerr,zRowH,zColH] = fastTxtRead(fname,indelim,headerRows,headerCols,skipRows,commentLine)
   delim = [];

   if (exist('indelim','var'))
       delim = regexprep(indelim,'\\t','\t');
   else
       delim = char(9);
   end

   if (exist('headerCols','var')==0 || isempty(headerCols))
       headerCols = 1;
   end

   if (exist('headerRows','var')==0 || isempty(headerRows))
       headerRows = 1;
   end

   if (exist('skipRows','var')==0 || isempty(skipRows))
       skipRows = 0;
   end

   if (exist('commentLine','var')==0 || isempty(commentLine))
    commentLine = [];
   end

    zdata = [];
    zerr = [];
    zRowH = [];
    zColH = [];


   if strcmp(computer,'MACI64')
       [zz,zout] = system([ '/usr/local/bin/greadlink -f ' fname]);
       zout = regexprep(zout,'\s*$','');
       if nargout > 2
         if isempty(commentLine)
            [zdata,zerr,zRowH,zColH] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols);
         else
            [zdata,zerr,zRowH,zColH] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols,commentLine);
         end
       else
         if isempty(commentLine)
            [zdata,zerr] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols);
         else
            [zdata,zerr] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols,commentLine);
         end
       end
   else
       [zz,zout] = system([ 'readlink -f ' zout]);
       zout = regexprep(zout,'\s*$','');
       if nargout > 2
         if isempty(commentLine)
            [zdata,zerr,zRowH,zColH] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols);
         else
            [zdata,zerr,zRowH,zColH] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols,commentLine);
         end
       else
         if isempty(commentLine)
            [zdata,zerr] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols);
         else
            [zdata,zerr] = fastTableReadV2(1,zout,[],delim,skipRows,headerRows,headerCols,commentLine);
         end
       end
   end
   % zdata = zdata';
end