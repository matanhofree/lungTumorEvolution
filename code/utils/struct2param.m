function outCell = structparam(strIn)
    
   fnames = fieldnames(strIn);
   cellval = struct2cell(strIn);
   
   ln = length(fnames)*2;
   outCell = cell(1,ln);
   outCell(1:2:end) = fnames;
   outCell(2:2:end) = cellval;
end