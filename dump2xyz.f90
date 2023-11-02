CHARACTER(255)   ::dummychar,inputfile,outputfile
CHARACTER(3)    ::atomtype
INTEGER         ::jj,kk,ntype,nnumber,natoms,nconfs,dummyint,ntemp,mol
REAL*8          ::xlo,xhi,ylo,yhi,zlo,zhi,x,y,z,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz
INTEGER,DIMENSION(:),ALLOCATABLE       ::nntype
REAL*8,DIMENSION(:),ALLOCATABLE         ::xx,yy,zz

CALL getarg(1, inputfile)
call getarg(2, outputfile)
nconfs=100
!READ(*,*)natoms
!READ(*,*)inputfile
!inputfile="/media/yuri/Seagate Backup Plus Drive/interface/medium-water/300/1/positions.lammpstrj"
!outputfile='/media/yuri/Seagate Backup Plus Drive/interface/medium-water/300/1/medium-water-300-1.xyz'
!READ(*,*)nconfs
!READ(*,*)ntemp

OPEN(9,FILE=trim(inputfile),STATUS='old')
READ(9,*)dummychar
READ(9,*)dummyint
READ(9,*)dummychar
READ(9,*)natoms
READ(9,*)dummychar
READ(9,*)xlo,xhi
READ(9,*)ylo,yhi
READ(9,*)zlo,zhi
READ(9,*)dummychar
do jj=1,natoms
  READ(9,*)nnumber,ntype,dummyint,x,y,z
enddo
close(9)
OPEN(9,FILE=trim(inputfile),STATUS='old')
OPEN(10,FILE=trim(outputfile))
DO kk=1,nconfs
  print*,kk
  !if ( mod(kk,100) == 0 ) print*,kk
  READ(9,*)dummychar
  READ(9,*)dummyint
  READ(9,*)dummychar
  READ(9,*)natoms
  READ(9,*)dummychar
  READ(9,*)xlo,xhi
  READ(9,*)ylo,yhi
  READ(9,*)zlo,zhi
  READ(9,*)dummychar
  WRITE(10,*)natoms
  WRITE(10,*)(xhi-xlo),(yhi-ylo),(zhi-zlo)
  DO jj=1,natoms
    READ(9,*)nnumber,ntype,dummyint,x,y,z
    SELECT CASE (ntype)
      CASE(1,4,6,8,10,12,23,31,32,37)
        atomtype='C'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (5,7,9,11,13,15,26,28,35,38)
        atomtype='H'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (2,14,19,25,27,29,30,33,34)
        atomtype='O'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (3)
        atomtype='N'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (16)
        atomtype='Cl'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (17)
        atomtype='Na'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (18,24,36)
        atomtype='S'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (20)
        atomtype='Mg'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (21)
        atomtype='Ca'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
      CASE (22)
        atomtype='K'
        WRITE(10,'(A3,3F9.4)')atomtype,x,y,z
    END SELECT
    !WRITE(10,'(A3,3F9.4)')atomtype,x,y,z  
  ENDDO
ENDDO
CLOSE(9)
CLOSE(10)

STOP
END
