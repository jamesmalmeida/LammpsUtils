program cuttrj

IMPLICIT NONE

double precision,allocatable :: coord(:,:)
double precision,dimension(3,2) :: box
double precision :: hcx,hcy,hcz,lcx,lcy,lcz
double precision :: radius,maxradius,minradius
integer :: i,j,k,l,nstep,natm,m,jump
integer :: l1,l2,l3
integer,allocatable :: idx(:),mol(:),type(:),iflag(:)
double precision,allocatable :: timeinregion(:,:),molcm(:,:)
character(LEN=30) :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dm(30),inputname
integer,allocatable :: mark(:,:),markpore(:,:,:)
integer,allocatable :: markporenumber(:)
double precision,dimension(2) :: minz,maxz
double precision :: distancex,distancey,distancez,binsize
integer :: nbins,timeinregionthr,filenum,maxmol
integer,allocatable :: maxbin(:),minbin(:),numatmbin(:),molatmnum(:)
character(LEN=3),allocatable :: atmlabel(:)
character(LEN=50),allocatable :: filename(:)

CALL getarg(1,inputname)
nstep=4000
nbins=80
timeinregionthr=0.5

open(20,file=inputname)

read(20,*)dummy1
read(20,*)dummy1
read(20,*)dummy1
read(20,*)natm
read(20,*)dummy1
read(20,*)box(1,1),box(1,2)
read(20,*)box(2,1),box(2,2)
read(20,*)box(3,1),box(3,2)
read(20,*)dummy1

allocate(coord(natm,3))
allocate(numatmbin(nbins),type(natm),filename(nbins+20))
allocate(markporenumber(nbins))
allocate(atmlabel(natm),mol(natm))

maxmol=0

do i=1,natm
  read(20,*)dm(1),type(i),mol(i),dm(2),dm(4),dm(5)
  if(mol(i).gt.maxmol)maxmol=mol(i)
  print*,i
enddo

allocate(molatmnum(maxmol),molcm(maxmol,3))
allocate(timeinregion(maxmol,nbins))
allocate(mark(maxmol,nstep))

rewind(20)

binsize=(box(3,2)-box(3,1))/nbins

write(*,*)natm,' atoms'
write(*,*)maxmol,' molecules'
write(*,*)box(1,1),box(1,2)
write(*,*)box(2,1),box(2,2)
write(*,*)box(3,1),box(3,2)

mark=1
markpore=1
markporenumber=0

l=0

do j=1,nstep
  molatmnum=0
  molcm=0
  print*,j
  read(20,*)dummy1
  read(20,*)dummy1
  read(20,*)dm(1)
  read(20,*)natm
  read(20,*)dummy1
  read(20,*)box(1,1),box(1,2)
  read(20,*)box(2,1),box(2,2)
  read(20,*)box(3,1),box(3,2)
  read(20,*)dm(1)
  m=0
  markporenumber=0
  do i=1,natm
    read(20,*)dm(1),dm(2),mol(i),(coord(i,k),k=1,3)
  enddo
!$OMP  PARALLEL DO &
!$OMP& DEFAULT(SHARED) PRIVATE(l,i,k) &
!$OMP& REDUCTION(+:MOLCM) &
!$OMP& REDUCTION(+:MOLATMNUM)
  do l=1,maxmol
    do i=1,natm
      if(mol(i).eq.l)then
        molcm(l,1)=molcm(l,1)+coord(i,1)
        molcm(l,2)=molcm(l,2)+coord(i,2)
        molcm(l,3)=molcm(l,3)+coord(i,3)
        molatmnum(l)=molatmnum(l)+1
      endif
    enddo 
    molcm(l,1)=molcm(l,1)/molatmnum(l)
    molcm(l,2)=molcm(l,2)/molatmnum(l)
    molcm(l,3)=molcm(l,3)/molatmnum(l)
    do k=1,nbins
      if(molcm(l,3).ge.((k-1)*binsize).and.molcm(l,3).lt.(k*binsize))then
        mark(l,j)=k
      endif        
    enddo
  enddo
!$OMP  END PARALLEL DO  
enddo

timeinregion=0

do j=1,nstep
  do l=1,maxmol
    do k=1,nbins
      if(mark(l,j)==k)timeinregion(l,k)=timeinregion(l,k)+1
    enddo 
  enddo
enddo

numatmbin=0

do l=1,maxmol
  do k=1,nbins
    timeinregion(l,k)=timeinregion(l,k)/nstep
    if(timeinregion(l,k)>timeinregionthr)then
      numatmbin(k)=numatmbin(k)+1*molatmnum(l)
    endif  
  enddo
enddo

do k=1,nbins
  filenum=20+k
  write(filename(filenum),'(A28,I4.4,A10)')"./regions-by-percentage/bin-",k,".lammpstrj"
  if(numatmbin(k)>0)open(filenum,file=filename(filenum))
enddo

rewind(20)

do j=1,nstep
  print*,'Writting Step ',j
  read(20,*)dm(1),dm(2)
  read(20,*)dm(3)
  read(20,*)dm(4),dm(5),dm(6),dm(7)
  read(20,*)natm
  read(20,*)dm(9),dm(10),dm(11),dm(12),dm(13),dm(14)
  read(20,*)box(1,1),box(1,2)
  read(20,*)box(2,1),box(2,2)
  read(20,*)box(3,1),box(3,2)
  read(20,*)dm(15),dm(16),dm(17),dm(18),dm(19),dm(20),dm(21),dm(22)
  do k=1,nbins
    if(numatmbin(k)>0)then
      write(20+k,'(A6,A9)')dm(1),dm(2)
      write(20+k,'(A8)')dm(3)
      write(20+k,'(A6,A7,A3,A6)')dm(4),dm(5),dm(6),dm(7)
      write(20+k,'(I6)')numatmbin(k)
      write(20+k,'(A6,A4,A7,3A3)')dm(9),dm(10),dm(11),dm(12),dm(13),dm(14)
      write(20+k,'(2F18.8)')box(1,1),box(1,2)
      write(20+k,'(2F18.8)')box(2,1),box(2,2)
      write(20+k,'(2F18.8)')box(3,1),box(3,2)
      write(20+k,'(A6,A6,A3,A5,A4,3A2)')dm(15),dm(16),dm(17),dm(18),dm(19),dm(20),dm(21),dm(22)
    endif
  enddo
  do i=1,natm
    read(20,*)dm(1),dm(2),mol(i),(coord(i,k),k=1,3)    
    do k=1,nbins
      if(numatmbin(k)>0)then
        if(timeinregion(mol(i),k)>timeinregionthr)then
          write(20+k,'(A7,A5,I8,3F18.8)')dm(1),dm(2),mol(i),coord(i,1),coord(i,2),coord(i,3)
        endif
      endif
    enddo
  enddo
enddo

end program cuttrj

