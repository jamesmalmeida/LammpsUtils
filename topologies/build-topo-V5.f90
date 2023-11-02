program buildtopo

IMPLICIT NONE

integer                       :: i,j,k,l,m
integer                       :: numberoffragments,temp,temp1,temp2,temp3,temp4
integer,allocatable           :: numberofatoms(:),numberofbonds(:),numberofangles(:),numberofdihedrals(:),numberofimpropers(:)
integer,allocatable           :: id(:,:),molecule(:,:),atomtype(:,:)
integer,allocatable           :: bondid(:,:),bondtype(:,:),bondatom(:,:,:)
integer,allocatable           :: angleid(:,:),angletype(:,:),angleatom(:,:,:)
integer,allocatable           :: dihedralid(:,:),dihedraltype(:,:),dihedralatom(:,:,:)
integer,allocatable           :: improperid(:,:),impropertype(:,:),improperatom(:,:,:)
double precision,allocatable  :: charge(:,:)
character(len=50),allocatable :: file(:)
character(len=50)             :: dummy,packfile

integer                       :: numberofatomspackmol,moleculetempcounter
integer                       :: totalnumberofatoms,totalnumberofbonds,totalnumberofangles
integer                       :: totalnumberofdihedrals,totalnumberofimpropers
integer,allocatable           :: numberofmolecules(:)
double precision              :: coordinates(3)

open(87,file='build-topo-V5.in')

read(87,*)numberoffragments

allocate(file(numberoffragments))
allocate(numberofatoms(numberoffragments))
allocate(numberofbonds(numberoffragments))
allocate(numberofangles(numberoffragments))
allocate(numberofdihedrals(numberoffragments))
allocate(numberofimpropers(numberoffragments))
allocate(numberofmolecules(numberoffragments))

do i=1,numberoffragments
  read(87,*)file(i)
enddo

read(87,*)packfile
open(30,file=packfile)

do i=1,numberoffragments
  read(87,*)numberofmolecules(i)
  print*,numberofmolecules(i),file(i),' molecules'
enddo

write(*,'(A14,2A20)')'Opening files '

do i=1,numberoffragments
  open(10+i,file=file(i))
enddo

do i=1,numberoffragments
  read(10+i,*)
  read(10+i,*)
enddo

!READ THE HEADER, MUST SPECIFY EVEN IF ZERO, FOR EXAMPLE: '0 dihedrals'.
do i=1,numberoffragments
  read(10+i,*)numberofatoms(i)
  write(*,'(A5,A20,A5,I5,A6)')'File ',file(i),' has ',numberofatoms(i),' atoms'
  read(10+i,*)numberofbonds(i)
  write(*,'(A5,A20,A5,I5,A6)')'File ',file(i),' has ',numberofbonds(i),' bonds'
  read(10+i,*)numberofangles(i)
  write(*,'(A5,A20,A5,I5,A7)')'File ',file(i),' has ',numberofangles(i),' angles' 
  read(10+i,*)numberofdihedrals(i)
  write(*,'(A5,A20,A5,I5,A10)')'File ',file(i),' has ',numberofdihedrals(i),' dihedrals'
  read(10+i,*)numberofimpropers(i)
  write(*,'(A5,A20,A5,I5,A10)')'File ',file(i),' has ',numberofimpropers(i),' impropers'
  read(10+i,*)
  read(10+i,*)
  read(10+i,*)
enddo

allocate(id(numberoffragments,maxval(numberofatoms)))
allocate(molecule(numberoffragments,maxval(numberofatoms)))
allocate(atomtype(numberoffragments,maxval(numberofatoms)))
allocate(charge(numberoffragments,maxval(numberofatoms)))

allocate(bondid(numberoffragments,maxval(numberofbonds)))
allocate(bondtype(numberoffragments,maxval(numberofbonds)))
allocate(bondatom(numberoffragments,maxval(numberofbonds),2))

allocate(angleid(numberoffragments,maxval(numberofangles)))
allocate(angletype(numberoffragments,maxval(numberofangles)))
allocate(angleatom(numberoffragments,maxval(numberofangles),3))

allocate(dihedralid(numberoffragments,maxval(numberofdihedrals)))
allocate(dihedraltype(numberoffragments,maxval(numberofdihedrals)))
allocate(dihedralatom(numberoffragments,maxval(numberofdihedrals),4))

allocate(improperid(numberoffragments,maxval(numberofimpropers)))
allocate(impropertype(numberoffragments,maxval(numberofimpropers)))
allocate(improperatom(numberoffragments,maxval(numberofimpropers),4))

!READ ATOM INFO, IGNORE THE BOX CROSSING FLAG.
!NOTE THAT THE TYPES NUMBER, BOX SIZE, AND PARAMETERS SHOULD NOT BE ON THE TOPOLOGY FILES. YOU SHOULD BUILD A HEADER TO ADD THAT TO
!THE OUTPUT FILE. 

do i=1,numberoffragments
  write(*,*)'Reading atom info of ',file(i)
  do j=1,numberofatoms(i)
    read(10+i,*)id(i,j),molecule(i,j),atomtype(i,j),charge(i,j) 
    !write(*,*)i,id(i,j),molecule(i,j),charge(i,j)
  enddo
enddo

do i=1,numberoffragments
  if(numberofbonds(i).ne.0)then
    read(10+i,*)
    read(10+i,*)
    read(10+i,*)
    write(*,*)'Reading bond info of ',file(i)
    do j=1,numberofbonds(i)
      read(10+i,*)bondid(i,j),bondtype(i,j),bondatom(i,j,1),bondatom(i,j,2)
      !write(*,*)i,bondid(i,j),bondtype(i,j),bondatom(i,j,1),bondatom(i,j,2)
    enddo
  endif
enddo

do i=1,numberoffragments
  if(numberofangles(i).ne.0)then
    read(10+i,*)
    read(10+i,*)
    read(10+i,*)
    write(*,*)'Reading angle info of ',file(i)
    do j=1,numberofangles(i)
      read(10+i,*)angleid(i,j),angletype(i,j),angleatom(i,j,1),angleatom(i,j,2),angleatom(i,j,3)
      !write(*,*)i,angleid(i,j),angletype(i,j),angleatom(i,j,1),angleatom(i,j,2),angleatom(i,j,3)
    enddo
  endif
enddo

do i=1,numberoffragments
  if(numberofdihedrals(i).ne.0)then
    read(10+i,*)
    read(10+i,*)
    read(10+i,*)          
    write(*,*)'Reading dihedral info of ',file(i)
    do j=1,numberofdihedrals(i)
      read(10+i,*)dihedralid(i,j),dihedraltype(i,j),dihedralatom(i,j,1),dihedralatom(i,j,2),dihedralatom(i,j,3),dihedralatom(i,j,4)
      !write(*,*)i,dihedralid(i,j),dihedraltype(i,j),dihedralatom(i,j,1),dihedralatom(i,j,2),dihedralatom(i,j,3),dihedralatom(i,j,4)
    enddo
  endif
enddo

do i=1,numberoffragments
  if(numberofimpropers(i).ne.0)then
    read(10+i,*)
    read(10+i,*)
    read(10+i,*)          
    write(*,*)'Reading improper info of ',file(i)
    do j=1,numberofimpropers(i)
      read(10+i,*)improperid(i,j),impropertype(i,j),improperatom(i,j,1),improperatom(i,j,2),improperatom(i,j,3),improperatom(i,j,4)
      !write(*,*)i,improperid(i,j),impropertype(i,j),improperatom(i,j,1),improperatom(i,j,2),improperatom(i,j,3),improperatom(i,j,4)
    enddo
  endif
enddo

print*,'Opening Packmol File'
read(30,*)numberofatomspackmol
read(30,*) !(numberofmolecules(k),k=1,numberoffragments)
print*,'Found ',numberofatomspackmol,' atoms on packmol file'
print*,'The number of molecules for each fragment is: ',numberofmolecules(:)

totalnumberofatoms=0
do i=1,numberoffragments
  totalnumberofatoms=totalnumberofatoms+numberofatoms(i)*numberofmolecules(i)
enddo

!print*,totalnumberofatoms

if(numberofatomspackmol.ne.totalnumberofatoms)then
  print*,'Error in the number of atoms or molecules in packmol xyz file'
  call exit
endif

totalnumberofbonds=0
do i=1,numberoffragments
  totalnumberofbonds=totalnumberofbonds+numberofbonds(i)*numberofmolecules(i)
enddo
totalnumberofangles=0
do i=1,numberoffragments
  totalnumberofangles=totalnumberofangles+numberofangles(i)*numberofmolecules(i)
enddo
totalnumberofdihedrals=0
do i=1,numberoffragments
  totalnumberofdihedrals=totalnumberofdihedrals+numberofdihedrals(i)*numberofmolecules(i)
enddo
totalnumberofimpropers=0
do i=1,numberoffragments
  totalnumberofimpropers=totalnumberofimpropers+numberofimpropers(i)*numberofmolecules(i)
enddo

open(40,file='HT.top')

write(40,*)
write(40,*)
write(40,*)totalnumberofatoms,' atoms'
write(40,*)totalnumberofbonds,' bonds'
write(40,*)totalnumberofangles,' angles'
write(40,*)totalnumberofdihedrals,' dihedrals'
write(40,*)totalnumberofimpropers,' impropers'
write(40,*)
write(40,*)'Atoms'
write(40,*)
moleculetempcounter=0
do i=1,numberoffragments
  do j=1,numberofmolecules(i)
    moleculetempcounter=moleculetempcounter+1
    do k=1,numberofatoms(i)
       read(30,*)dummy,(coordinates(l),l=1,3)
       temp=id(i,k)+numberofatoms(i)*(j-1)
       do m=1,i-1
         temp=temp+numberofatoms(i-m)*numberofmolecules(i-m)
       enddo
!       write(40,'(I7,I6,I4,F10.4,3F15.8)') temp,molecule(i,k),atomtype(i,k),charge(i,k),coordinates(:)
       write(40,'(I7,I6,I4,F10.4,3F15.8)') temp,moleculetempcounter,atomtype(i,k),charge(i,k),coordinates(:)
    enddo
  enddo
enddo
write(40,*)
write(40,*)'Bonds'
write(40,*)
do i=1,numberoffragments
  if(numberofbonds(i).ne.0)then
    do j=1,numberofmolecules(i)
      do k=1,numberofbonds(i)
      temp=bondid(i,k)+numberofbonds(i)*(j-1)
      temp1=bondatom(i,k,1)+numberofatoms(i)*(j-1)
      temp2=bondatom(i,k,2)+numberofatoms(i)*(j-1)
      do m=1,i-1
        temp=temp+numberofbonds(i-m)*numberofmolecules(i-m)
        temp1=temp1+numberofatoms(i-m)*numberofmolecules(i-m)
        temp2=temp2+numberofatoms(i-m)*numberofmolecules(i-m)
      enddo
      write(40,*)temp,bondtype(i,k),temp1,temp2
      enddo
    enddo
  endif
enddo
write(40,*)
write(40,*)'Angles'
write(40,*)
do i=1,numberoffragments
  if(numberofangles(i).ne.0)then
    do j=1,numberofmolecules(i)
      do k=1,numberofangles(i)
      temp=angleid(i,k)+numberofangles(i)*(j-1)
      temp1=angleatom(i,k,1)+numberofatoms(i)*(j-1)
      temp2=angleatom(i,k,2)+numberofatoms(i)*(j-1)
      temp3=angleatom(i,k,3)+numberofatoms(i)*(j-1)
      do m=1,i-1
        temp=temp+numberofangles(i-m)*numberofmolecules(i-m)
        temp1=temp1+numberofatoms(i-m)*numberofmolecules(i-m)
        temp2=temp2+numberofatoms(i-m)*numberofmolecules(i-m)
        temp3=temp3+numberofatoms(i-m)*numberofmolecules(i-m)
      enddo
      write(40,*)temp,angletype(i,k),temp1,temp2,temp3
      enddo
    enddo
  endif
enddo
write(40,*)
write(40,*)'Dihedrals'
write(40,*)
do i=1,numberoffragments
  if(numberofdihedrals(i).ne.0)then
    do j=1,numberofmolecules(i)
      do k=1,numberofdihedrals(i)
        temp=dihedralid(i,k)+numberofdihedrals(i)*(j-1)
        temp1=dihedralatom(i,k,1)+numberofatoms(i)*(j-1)
        temp2=dihedralatom(i,k,2)+numberofatoms(i)*(j-1)
        temp3=dihedralatom(i,k,3)+numberofatoms(i)*(j-1)
        temp4=dihedralatom(i,k,4)+numberofatoms(i)*(j-1)
        do m=1,i-1
          temp=temp+numberofdihedrals(i-m)*numberofmolecules(i-m)
          temp1=temp1+numberofatoms(i-m)*numberofmolecules(i-m)
          temp2=temp2+numberofatoms(i-m)*numberofmolecules(i-m)
          temp3=temp3+numberofatoms(i-m)*numberofmolecules(i-m)
          temp4=temp4+numberofatoms(i-m)*numberofmolecules(i-m)
        enddo
        write(40,*)temp,dihedraltype(i,k),temp1,temp2,temp3,temp4
      enddo
    enddo
  endif
enddo
write(40,*)
write(40,*)'Impropers'
write(40,*)
do i=1,numberoffragments
  if(numberofimpropers(i).ne.0)then
    do j=1,numberofmolecules(i)
      do k=1,numberofimpropers(i)
        temp=improperid(i,k)+numberofimpropers(i)*(j-1)
        temp1=improperatom(i,k,1)+numberofatoms(i)*(j-1)
        temp2=improperatom(i,k,2)+numberofatoms(i)*(j-1)
        temp3=improperatom(i,k,3)+numberofatoms(i)*(j-1)
        temp4=improperatom(i,k,4)+numberofatoms(i)*(j-1)
        do m=1,i-1
          temp=temp+numberofimpropers(i-m)*numberofmolecules(i-m)
          temp1=temp1+numberofatoms(i-m)*numberofmolecules(i-m)
          temp2=temp2+numberofatoms(i-m)*numberofmolecules(i-m)
          temp3=temp3+numberofatoms(i-m)*numberofmolecules(i-m)
          temp4=temp4+numberofatoms(i-m)*numberofmolecules(i-m)
        enddo
        write(40,*)temp,impropertype(i,k),temp1,temp2,temp3,temp4
      enddo
    enddo
  endif
enddo

end program buildtopo
