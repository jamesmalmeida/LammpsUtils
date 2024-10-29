program identify

implicit none

integer                      :: i,j,k,natoms
integer,allocatable          :: bondcounter(:),bondindex(:,:)
double precision,allocatable :: coordinates(:,:),charge(:)
double precision             :: bondcutoff,distance
character,allocatable        :: atomsymbol(:),bondtype(:,:)
character(len=5),allocatable :: bonds(:),charmmtype(:)

open(10,file='gaussian.xyz')

bondcutoff=1.7

read(10,*)natoms
read(10,*)

allocate(coordinates(natoms,3),atomsymbol(natoms))
allocate(bondcounter(natoms),bondtype(natoms,5),bondindex(natoms,5))
allocate(charmmtype(natoms),bonds(natoms),charge(natoms))


do i=1,natoms
  read(10,*)atomsymbol(i),(coordinates(i,k),k=1,3),charge(i)
enddo

bondcounter=0

do i=1,natoms
  do j=i+1,natoms
    distance=sqrt(((coordinates(i,1)-coordinates(j,1))**2)+&
                 &((coordinates(i,2)-coordinates(j,2))**2)+&
                 &((coordinates(i,3)-coordinates(j,3))**2))
    if(distance.le.bondcutoff)then
      bondcounter(i)=bondcounter(i)+1
      bondcounter(j)=bondcounter(j)+1
      bondtype(i,bondcounter(i))=atomsymbol(j)
      bondtype(j,bondcounter(j))=atomsymbol(i)
      bondindex(i,bondcounter(i))=j
      bondindex(j,bondcounter(j))=i
    endif
  enddo
enddo

do i=1,natoms
  write(bonds(i),'(5A1)')atomsymbol(i),(bondtype(i,j),j=1,bondcounter(i))
  write(*,'(A5,I4,A4,A5)')'Atom ',i,' is ',bonds(i)
enddo
do i=1,natoms
  if(bonds(i)=="HC   ")then
    if(bonds(bondindex(i,bondcounter(i)))=="CCCHH ")charmmtype(i)="HGA2"
    if(bonds(bondindex(i,bondcounter(i)))=="CCHHC ")charmmtype(i)="HGA2"
    if(bonds(bondindex(i,bondcounter(i)))=="CHHCC ")charmmtype(i)="HGA2"
    if(bonds(bondindex(i,bondcounter(i)))=="CCHHH ")charmmtype(i)="HGA3"
    if(bonds(bondindex(i,bondcounter(i)))=="CHHHC ")charmmtype(i)="HGA3"
    if(bonds(bondindex(i,bondcounter(i)))=="CHHHH ")charmmtype(i)="HGA3"
  endif
  SELECT CASE (bonds(i))
     CASE ("OCH  ")
        charmmtype(i)='OG311'
     CASE ("OC   ")
        charmmtype(i)='OG2D1'
     CASE ("CCCHH","CCHHC","CHHCC")
        charmmtype(i)='CG321'
     CASE ("CCHHH","CHHHC","CHHHH")
        charmmtype(i)='CG331'
     CASE ("COOC ")
        charmmtype(i)='CG2O2'
     CASE ("HO   ")
        charmmtype(i)='HGP1'
     CASE("HC    ")

     CASE DEFAULT
        print*,'Problem With Atom',i
  END SELECT  
enddo

open(20,file='charmm.xyz')

write(20,*)natoms
write(20,*)'XYZ with Charmm atom types'
do i=1,natoms
  write(20,'(A6,4F15.8)')charmmtype(i),(coordinates(i,k),k=1,3),charge(i)
enddo

end program identify
