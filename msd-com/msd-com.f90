program msdcom
IMPLICIT NONE

double precision,allocatable :: coord(:,:),coord0(:,:),msd(:,:)
double precision             :: zlo,zhi,binsize
character(LEN=50)            :: inputfile,inputtraj,outputfile,dummy
integer                      :: i,j,k,l,nsteps,natoms,nbins,windowsize,nwindows
integer,allocatable          :: counter(:,:)

CALL getarg(1, inputfile)

OPEN(9,file=inputfile)

READ(9,*)inputtraj
READ(9,*)outputfile
READ(9,*)nbins
READ(9,*)nsteps
READ(9,*)windowsize
CLOSE(9)

OPEN(10,file=inputtraj)

zhi=-1000
zlo= 1000
READ(10,*)dummy
READ(10,*)dummy
READ(10,*)dummy
READ(10,*)dummy,natoms
allocate(coord(3,natoms),coord0(3,natoms),msd(nbins,windowsize),counter(nbins,windowsize))
do i=1,natoms
  READ(10,*)dummy,coord0(1,i),coord0(2,i),coord0(3,i)
  if(coord0(3,i)>zhi)zhi=coord0(3,i)
  if(coord0(3,i)<zlo)zlo=coord0(3,i)
enddo
close(10)

nwindows=nsteps/windowsize
binsize=(zhi-zlo)/nbins
print*,(zhi-zlo),binsize,nwindows

msd=0
counter=0
OPEN(10,file=inputtraj)
READ(10,*)dummy
READ(10,*)dummy
READ(10,*)dummy

do l=1,nwindows
  do j=1,windowsize
    print*,'Window ',l,'Step ',j
    read(10,*)dummy
    if(j==1)then
      do i=1,natoms
        READ(10,*)dummy,coord0(1,i),coord0(2,i),coord0(3,i)
        coord(1,i)=coord0(1,i)
        coord(2,i)=coord0(2,i)
        coord(3,i)=coord0(3,i)
      enddo
    else 
      do i=1,natoms
        READ(10,*)dummy,coord(1,i),coord(2,i),coord(3,i)
      enddo
    endif
    do i=1,natoms
      do k=1,nbins
        if(coord(3,i)>=(binsize*(k-1)).and.coord(3,i)<(binsize*k))then
          msd(k,j)=msd(k,j)+((coord(1,i)-coord0(1,i))**2+(coord(2,i)-coord0(2,i))**2+(coord(3,i)-coord0(3,i))**2)
          counter(k,j)=counter(k,j)+1
        endif
      enddo 
    enddo
  enddo
enddo
close(10)

open(20,file=outputfile)
do j=1,windowsize
  do k=1,nbins
    if(counter(k,j)<10)then
      do l=1,windowsize
        msd(k,l)=0
      enddo
    endif
    if(counter(k,j).ne.0)msd(k,j)=msd(k,j)/counter(k,j)
  enddo
  write(20,'(101F16.8)')(j*0.5*0.1),msd(:,j)
enddo
close(20)

deallocate(coord,coord0,msd,counter)

end program msdcom
