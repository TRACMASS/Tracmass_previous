
load decaymask;

aid=fopen('totnumtrajmatr.bin','a');

if (exist('rd.mat','file'))
  load rd;
end

%if (exist('totnumtraj.mat','file'))
%  load totnumtraj;
%end

mysql ('open','localhost','root','','traj')       

x=nc_varget('/data/GOM/grid.cdf','x');
dx=x(2:80,2:164)-x(1:79,1:163);
dx(80,:)=dx(79,:);

y=nc_varget('/data/GOM/grid.cdf','y');
dy=y(2:80,2:164)-y(1:79,1:163);
dy(80,:)=dy(79,:);

depth_s=-nc_varget('/data/GOM/fields/grid.cdf','depth');
depth_s=-nc_varget('/data/GOM/fields/grid.cdf','depth_at_sigma');

maxdepth_s= max(max(max(depth_s)));

minIntsStruct=mysql('SELECT min(ints) FROM trajgoKB01');
minInts=minIntsStruct.min_ints_;

maxIntsStruct=mysql('SELECT max(ints) FROM trajgoKB01');
maxInts=maxIntsStruct.max_ints_;


for    ttt=1:300
  
  trajres=mysql(['SELECT count(round(x1)) countraj,round(x1) x1, round(y1)  y1' ...
                 ' FROM trajgoKB01 WHERE  ints=' ...
                 num2str(minInts+ttt-1) ...
                 ' AND x1>1 AND y1>1 AND x1<163 AND y1<79 ' ...
                 'GROUP BY round(x1), round(y1);']); 
  
  %disp([ttt length(trajres)])
  
  x1=0; y1=0;
  y1=[trajres.x1];
  x1=[trajres.y1];  
  countraj=[trajres.countraj]; 
    
  numtraj=x.*NaN;
  for ij=1:length(x1)
    numtraj(x1(ij),y1(ij))=countraj(ij);
  end
  resdecay(ttt)=nansum(nansum(numtraj.*squeeze(decaymask(1,:,:))));
%  timenumtraj(ttt,:,:)=numtraj;
fwrite (aid,numtraj,'single');
end

if (~exist('rd','var'))
  rd=resdecay;
else
  rd(size(rd,1)+1,:)=resdecay;
end
save rd rd

%if (~exist('totnumtraj','var'))
%  totnumtraj(1,:,:,:)=timenumtraj;
%else
 % totnumtraj(size(totnumtraj,1),:,:,:)=timenumtraj;
%end
%save totnumtraj totnumtraj

quit