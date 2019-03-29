function DBM=FDK(nw,TFL)

%... Construction of Drained Blocks
DBM=zeros(size(TFL));

index=zeros(numel(TFL),1);
TFLV=sort(TFL(:));
WellLocInd=find(~TFL);
[nx,ny, nz]=size(TFL);

for i=1:nw;
  DBM(WellLocInd(i))=i;  
  index(i)=WellLocInd(i);
end
for k=nw+1:numel(TFL)
%   clear index 
  index=find(TFL==TFLV(k));
  
 for ttt=1:numel(index)
  %NBELIND=zeros(4,1);
  rgtelind=index(ttt)+nx;   %%% rightside element index
  lftelind=index(ttt)-nx;   %%% leftside element index
  upelind=index(ttt)-1;     %%% upside element index
  dwnelind=index(ttt)+1;    %%% downside element index
  
  if ~isempty(intersect(dwnelind,1:nx:nx*ny))
      dwnelind=NaN;
  elseif ~isempty(intersect(upelind,nx:nx:nx*ny))
       upelind=NaN;
  end
  
  NBELIND=[rgtelind; lftelind; upelind; dwnelind];
  NBELINDind= NBELIND>0 & NBELIND<nx*ny+1;

  NGBL=unique(nonzeros(DBM(NBELIND(NBELINDind))),'rows');
 
  NGBL=NGBL(NGBL~=8);
%DBM(index(ttt))=NGBL(1);

   if size(NGBL,1)==1
      DBM(index(ttt))=NGBL;
   else
       DBM(index(ttt))=8;
   end

end
%   pause(0.0001)
%  imagesc(DBM)
end

end