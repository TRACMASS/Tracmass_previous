
program extract
  
  PARAMETER (IM=164,JM=80,KB=22)
  
  real dx(im,jm), dy(im,jm), dz(kb)
  
open (41,file='/data/GOM/fort.41')

  do j=1,jm
     read(41,'(f10.3)') (dx(i,j),dy(i,j),i=1,im,1)
  enddo
  
  read(41,'(f10.3)') dz
  
print *,dz

end program extract
