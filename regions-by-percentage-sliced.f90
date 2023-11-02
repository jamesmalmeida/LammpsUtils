program cuttrj

double precision,allocatable :: coord(:,:,:)
double precision,dimension(3,2) :: box
double precision :: hcx,hcy,hcz,lcx,lcy,lcz
double precision :: radius,maxradius,minradius
integer :: i,j,k,l,nstep,natm,m
integer,allocatable :: idx(:),mol(:),type(:,:),iflag(:)
double precision,allocatable :: timeinregion(:,:)
character(LEN=30) :: dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dm(9)
integer,allocatable :: mark(:,:),markpore(:,:,:)
integer,dimension(3) :: markporenumber
double precision,dimension(2) :: minz,maxz
double precision :: distancex,distancey,distancez
double precision :: threashold
integer :: bla1,bla2,centerofmasscount,slicesize
character(len=50) :: filename(6)
double precision :: centerofmass(2)

threashold=0.9
slicesize=1
nstep=800

open(20,file='dump-production.lammpstrj')

read(20,*)dummy1
read(20,*)dummy1
read(20,*)dummy1
read(20,*)natm
read(20,*)dummy1
read(20,*)box(1,1),box(1,2)
read(20,*)box(2,1),box(2,2)
read(20,*)box(3,1),box(3,2)
read(20,*)dummy1
rewind 20

allocate(coord(1:natm,9,1:nstep),idx(1:natm),mol(1:natm),iflag(1:natm),type(1:natm,1:nstep))
allocate(markpore(1:natm,1:nstep,1:3),mark(1:natm,1:nstep))
allocate(timeinregion(1:natm,1:3))

do j=1,nstep
  print*,j
!  centerofmass=0
!  centerofmasscount=0
  read(20,*)dummy1
  read(20,*)dummy1
  read(20,*)dummy1
  read(20,*)dummy1
  read(20,*)dummy1
  read(20,*)box(1,1),box(1,2)
  read(20,*)box(2,1),box(2,2)
  read(20,*)box(3,1),box(3,2)
  read(20,*)dummy1
  do i=1,natm
    read(20,*)idx(i),type(i,j),dummy1,(coord(i,k,j),k=1,3)
!    if(coord(i,1,j)>((box(1,2)-box(1,1))/2))coord(i,1,j)=coord(i,1,j)-(box(1,2)-box(1,1))
!    if(coord(i,2,j)>((box(2,2)-box(2,1))/2))coord(i,2,j)=coord(i,2,j)-(box(2,2)-box(2,1))
!    if(coord(i,3,j)<=(30).and.coord(i,3,j)>=(-30).and.type(i,j)==7)then
!      centerofmass(1)=centerofmass(1)+coord(i,1,j)
!      centerofmass(2)=centerofmass(2)+coord(i,2,j)
!      centerofmasscount=centerofmasscount+1
!    endif
  enddo
!  centerofmass(1)=centerofmass(1)/centerofmasscount
!  centerofmass(2)=centerofmass(2)/centerofmasscount
!  write(*,'(2F12.6)'),centerofmass(:)
enddo

m=0
markporenumber=0
hcz=30
lcz=-30

do bla1=1,nstep-slicesize,slicesize
  markporenumber=0
  markpore=0
  mark=0
  maxradius=-1000
  minradius=1000

  bla2=bla1+slicesize
  print*,bla1,bla2

!  do i=1,natm
!    if(coord(i,3,bla1)<=(30).and.coord(i,3,bla1)>=(-30).and.type(i,bla1)==1)then
!      radius=sqrt(coord(i,1,1)**2+coord(i,2,1)**2)
!      if(radius>=maxradius)maxradius=radius
!      if(radius<=minradius)minradius=radius
!    endif
!  enddo  

  do j=bla1,bla2
    do i=1,natm
      if((coord(i,3,j)<=(hcz).and.coord(i,3,j)>=(lcz)))then
          radius=(coord(i,3,j)
          if(radius<=(minradius))then
            markpore(i,j,1)=1
            markpore(i,j,2)=0
            markpore(i,j,3)=0
            markporenumber(1)=markporenumber(1)+1
            elseif(radius<=(minradius).and.radius>=(minradius-5))then
            markpore(i,j,1)=0
            markpore(i,j,2)=1
            markpore(i,j,3)=0
            markporenumber(2)=markporenumber(2)+1
            elseif(radius>=(minradius))then
            markpore(i,j,1)=0
            markpore(i,j,2)=0
            markpore(i,j,3)=1
            markporenumber(3)=markporenumber(3)+1
          else
            markpore(i,j,:)=0
          endif
        else
          markpore(i,j,:)=0
        endif
        m=m+1
      endif
    enddo
  enddo
 
  print*,markporenumber(:)/(bla2-bla1)
 
  timeinregion=0
  do j=bla1,bla2
    do i=1,natm
      if(markpore(i,j,1)==1)timeinregion(i,1)=timeinregion(i,1)+1
      if(markpore(i,j,2)==1)timeinregion(i,2)=timeinregion(i,2)+1
      if(markpore(i,j,3)==1)timeinregion(i,3)=timeinregion(i,3)+1
    enddo
  enddo
  
  do i=1,natm
    if(type(i,j)>=7)then
      timeinregion(i,1)=timeinregion(i,1)/(bla2-bla1) 
      timeinregion(i,2)=timeinregion(i,2)/(bla2-bla1)
      timeinregion(i,3)=timeinregion(i,3)/(bla2-bla1)
    endif
  enddo
  
  !DIVIDE TIMEINREGION PER NUMBER OF STEPS
  !PRINT PER ATOM TIME IN REGION TO CHECK
  !SET THREASHOLD OF PERMANENCY IN A REGION TO CONSIDER THE ATOM BELONGING THERE
  !PRINT CORRESPONDING TO THAT
  
  m=0
  markporenumber=0
  do i=1,natm
    if(timeinregion(i,1)>=threashold)markporenumber(1)=markporenumber(1)+1
    if(timeinregion(i,2)>=threashold)markporenumber(2)=markporenumber(2)+1
    if(timeinregion(i,3)>=threashold)markporenumber(3)=markporenumber(3)+1
  enddo
  
  print*,'FINAL NUMBER OF ATOMS SHELL1',markporenumber(1)
  print*,'FINAL NUMBER OF ATOMS SHELL2',markporenumber(2)
  print*,'FINAL NUMBER OF ATOMS SHELL3',markporenumber(3)

  write(filename(1),'(A18,I4.4,A4)')"./regions/region1-",bla2,".xyz"
  write(filename(2),'(A18,I4.4,A4)')"./regions/region2-",bla2,".xyz"
  write(filename(3),'(A18,I4.4,A4)')"./regions/region3-",bla2,".xyz"
  write(filename(4),'(A21,I4.4,A4)')"./regions/region1-gr-",bla2,".xyz"
  write(filename(5),'(A21,I4.4,A4)')"./regions/region2-gr-",bla2,".xyz"
  write(filename(6),'(A21,I4.4,A4)')"./regions/region3-gr-",bla2,".xyz"
  !open(110,file=filename(1))
  !open(120,file=filename(2))
  !open(130,file=filename(3))
  open(140,file=filename(4))
  open(150,file=filename(5))
  open(160,file=filename(6))
  
  do j=bla1,bla2
     !write(110,*)markporenumber(1)
     !write(110,*)
     write(140,*)markporenumber(1)
     write(140,*)
  
     !write(120,*)markporenumber(2)
     !write(120,*)
     write(150,*)markporenumber(2)
     write(150,*)
  
     !write(130,*)markporenumber(3)
     !write(130,*)
     write(160,*)markporenumber(3)
     write(160,*)
  
     do i=1,natm
       if(timeinregion(i,1)>=threashold)then
         !select case(type(i,j))
         ! case(7)
         !  write(110,'(A3,2A10,F12.6)')'C  ','0.0000000','0.0000000',coord(i,3,j)
         ! case(8)
         !  write(110,'(A3,2A10,F12.6)')'O  ','0.0000000','0.0000000',coord(i,3,j)
         ! case(9)
         !  write(110,'(A3,2A10,F12.6)')'Na ','0.0000000','0.0000000',coord(i,3,j)
         ! case(10,12)
         !  write(110,'(A3,2A10,F12.6)')'Cl ','0.0000000','0.0000000',coord(i,3,j)
         ! case(11)
         !  write(110,'(A3,2A10,F12.6)')'Ca ','0.0000000','0.0000000',coord(i,3,j)
         !end select
         select case(type(i,j))
          case(7)
           write(140,'(A3,3F16.6)')'C  ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(8)
           write(140,'(A3,3F16.6)')'O  ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(9)
           write(140,'(A3,3F16.6)')'Na ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(10,12)
           write(140,'(A3,3F16.6)')'Cl ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(11)
           write(140,'(A3,3F16.6)')'Ca ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
         end select
       endif
       if(timeinregion(i,2)>=threashold)then
         !select case(type(i,j))
         ! case(7)
         !  write(120,'(A3,2A10,F12.6)')'C  ','0.0000000','0.0000000',coord(i,3,j)
         ! case(8)
         !  write(120,'(A3,2A10,F12.6)')'O  ','0.0000000','0.0000000',coord(i,3,j)
         ! case(9)
         !  write(120,'(A3,2A10,F12.6)')'Na ','0.0000000','0.0000000',coord(i,3,j)
         ! case(10,12)
         !  write(120,'(A3,2A10,F12.6)')'Cl ','0.0000000','0.0000000',coord(i,3,j)
         ! case(11)
         !  write(120,'(A3,2A10,F12.6)')'Ca ','0.0000000','0.0000000',coord(i,3,j)
         !end select
         select case(type(i,j))
          case(7)
           write(150,'(A3,3F16.6)')'C  ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(8)
           write(150,'(A3,3F16.6)')'O  ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(9)
           write(150,'(A3,3F16.6)')'Na ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(10,12)
           write(150,'(A3,3F16.6)')'Cl ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(11)
           write(150,'(A3,3F16.6)')'Ca ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
         end select
       endif
       if(timeinregion(i,3)>=threashold)then
         !select case(type(i,j))
         ! case(7)
         !  write(130,'(A3,2A10,F12.6)')'C  ','0.0000000','0.0000000',coord(i,3,j)
         ! case(8)
         !  write(130,'(A3,2A10,F12.6)')'O  ','0.0000000','0.0000000',coord(i,3,j)
         ! case(9)
         !  write(130,'(A3,2A10,F12.6)')'Na ','0.0000000','0.0000000',coord(i,3,j)
         ! case(10,12)
         !  write(130,'(A3,2A10,F12.6)')'Cl ','0.0000000','0.0000000',coord(i,3,j)
         ! case(11)
         !  write(130,'(A3,2A10,F12.6)')'Ca ','0.0000000','0.0000000',coord(i,3,j)
         !end select
         select case(type(i,j))
          case(7)
           write(160,'(A3,3F16.6)')'C  ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(8)
           write(160,'(A3,3F16.6)')'O  ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(9)
           write(160,'(A3,3F16.6)')'Na ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(10,12)
           write(160,'(A3,3F16.6)')'Cl ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
          case(11)
           write(160,'(A3,3F16.6)')'Ca ',coord(i,1,j),coord(i,2,j),coord(i,3,j)
         end select
       endif
     enddo
  enddo

  !close(110)
  !close(120)
  !close(130)
  close(140)
  close(150)
  close(160)
enddo

end program cuttrj

