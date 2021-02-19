function [zABC,ia,ib,ic] = intersect3(setA,setB,setC)
    [zAB] = intersect(setA,setB);
    [zABC] = intersect(zAB,setC);
    
    [zz,~,ia] = intersect(zABC,setA);
    assert(isequal(zz,zABC));
    [zz,~,ib] = intersect(zABC,setB);
    assert(isequal(zz,zABC));
    [zz,~,ic] = intersect(zABC,setC);
    assert(isequal(zz,zABC));
end