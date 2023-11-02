program hbonds

use omp_lib

IMPLICIT NONE

integer                      :: i,j,k,l,m,natoms,wateroxygentype,waterhydrogentype,h1,h2,h2onum
integer                      :: nstep,hbnum,hmol
integer,allocatable          :: id(:),type(:),mol(:)
double precision,allocatable :: pos(:,:) 
double precision             :: lox,loy,loz,hix,hiy,hiz,dist,dist1,dist2,theta,thetacut
double precision             :: oodistcut,v1(3),v2(3),v1m,v2m,v1v2
double precision             :: binsize
integer                      :: nbins
integer,allocatable          :: h2obinnum(:),hbbinnum(:)
integer                      :: nthreads,threadid
double precision             :: charge,start,finish,t1,delta
double precision             :: distx,disty,distz,dist1x,dist1y,dist1z,dist2x,dist2y,dist2z

nstep=700
oodistcut=4.0
wateroxygentype=14
waterhydrogentype=15
thetacut=120
thetacut=thetacut*0.0174533
binsize=1

open(9,file='positions.lammpstrj')

h2onum=0

read(9,*)
read(9,*)
read(9,*)
read(9,*)natoms
read(9,*)
read(9,*)lox,hix
read(9,*)loy,hiy
read(9,*)loz,hiz
read(9,*)
allocate(pos(3,natoms),id(natoms),type(natoms),mol(natoms))
do i=1,natoms
  read(9,*)id(j),type(j),mol(j),charge,pos(1,j),pos(2,j),pos(3,j)
  if(type(j)==14)h2onum=h2onum+1
enddo
rewind(9)

nbins=1+(hiz-loz)/binsize
allocate(h2obinnum(nbins),hbbinnum(nbins))

print*,h2onum,' Water Molecules'
print*,nbins,' bins'
h2obinnum=0
hbbinnum=0

do i=1,nstep
  call cpu_time(start)
  t1=secnds(0.0)
  print*,i
  hbnum=0
!  h2obinnum=0
!  hbbinnum=0
  read(9,*)
  read(9,*)
  read(9,*)
  read(9,*)natoms
  read(9,*)
  read(9,*)lox,hix
  read(9,*)loy,hiy
  read(9,*)loz,hiz
  read(9,*)
  do j=1,natoms
    read(9,*)id(j),type(j),mol(j),charge,pos(1,j),pos(2,j),pos(3,j)
    pos(3,j)=pos(3,j)-loz
  enddo  
  !$OMP  PARALLEL DO &
  !$OMP& default(none) &
  !$OMP& REDUCTION(+:hbbinnum) &
  !$OMP& REDUCTION(+:h2obinnum) &
  !$OMP& REDUCTION(+:hbnum) &
  !$OMP& REDUCTION(+:hmol) &
  !$OMP& PRIVATE(j,k,l,m) &
  !$OMP& PRIVATE(dist,dist1,dist2,h1,h2,v1,v2,v1m,v2m,v1v2,theta,threadid) &
  !$OMP& PRIVATE(distx,disty,distz,dist1x,dist1y,dist1z,dist2x,dist2y,dist2z) &
  !$OMP& SHARED(i,natoms,wateroxygentype,waterhydrogentype,nbins,pos,thetacut) &
  !$OMP& SHARED(mol,type,oodistcut,nthreads,hix,lox,hiy,loy,hiz,loz)
  do j=1,natoms
!    nthreads = OMP_GET_NUM_THREADS()
!    threadid = OMP_GET_THREAD_NUM()
!    print*,i,j,threadid,nthreads
    if(type(j)==wateroxygentype)then
      do k=1,nbins
        nthreads = OMP_GET_NUM_THREADS()
        threadid = OMP_GET_THREAD_NUM()
        !print*,k,h2obinnum(k),threadid,nthreads
        if(pos(3,j)>=(k-1).and.pos(3,j)<k)h2obinnum(k)=h2obinnum(k)+1
      enddo
      do k=1,natoms
        if(type(k)==wateroxygentype.and.j.ne.k)then
          distx=abs(pos(1,j)-pos(1,k))
          if(distx>((hix-lox)/2))distx=distx-(hix-lox)
          disty=abs(pos(2,j)-pos(2,k))
          if(disty>((hiy-loy)/2))disty=disty-(hiy-loy)
          distz=abs(pos(3,j)-pos(3,k))
          if(distz>((hiz-loz)/2))distz=distz-(hiz-loz)
          dist=sqrt(distx**2+disty**2+distz**2)
          if(dist.le.oodistcut)then
            hmol=0
            do l=1,natoms
              if(mol(k)==mol(l).and.type(l)==waterhydrogentype)then
                hmol=hmol+1
                if(hmol==1)then
                  h1=l
                  dist1x=abs(pos(1,j)-pos(1,l))
                  if(dist1x>((hix-lox)/2))dist1x=dist1x-(hix-lox)
                  dist1y=abs(pos(2,j)-pos(2,l))
                  if(dist1y>((hiy-loy)/2))dist1y=dist1y-(hiy-loy)
                  dist1z=abs(pos(3,j)-pos(3,l))
                  if(dist1z>((hiz-loz)/2))dist1z=dist1z-(hiz-loz)
                  dist1=sqrt(dist1x**2+dist1y**2+dist1z**2)
                endif
                if(hmol==2)then
                  h2=l
                  dist2x=abs(pos(1,j)-pos(1,l))
                  if(dist2x>((hix-lox)/2))dist2x=dist2x-(hix-lox)
                  dist2y=abs(pos(2,j)-pos(2,l))
                  if(dist2y>((hiy-loy)/2))dist2y=dist2y-(hiy-loy)
                  dist2z=abs(pos(3,j)-pos(3,l))
                  if(dist2z>((hiz-loz)/2))dist2z=dist2z-(hiz-loz)
                  dist2=sqrt(dist2x**2+dist2y**2+dist2z**2)
                  if(dist1.lt.dist2)then
                    v1(1)=pos(1,j)-pos(1,h1)
                    if(v1(1)>((hix-lox)/2))v1(1)=v1(1)-(hix-lox)
                    v1(2)=pos(2,j)-pos(2,h1)
                    if(v1(2)>((hiy-loy)/2))v1(2)=v1(2)-(hiy-loy)
                    v1(3)=pos(3,j)-pos(3,h1)
                    if(v1(3)>((hiz-loz)/2))v1(3)=v1(3)-(hiz-loz)
                    v2(1)=pos(1,k)-pos(1,h1)
                    if(v2(1)>((hix-lox)/2))v2(1)=v2(1)-(hix-lox)
                    v2(2)=pos(2,k)-pos(2,h1)
                    if(v2(2)>((hiy-loy)/2))v2(2)=v2(2)-(hiy-loy)
                    v2(3)=pos(3,k)-pos(3,h1)
                    if(v2(3)>((hiz-loz)/2))v2(3)=v2(3)-(hiz-loz)
                    v1m=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
                    v2m=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
                    v1v2=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
                    theta=acos( v1v2/(v1m*v2m) )
                    if(theta.ge.thetacut)then
                      do m=1,nbins
                        if(pos(3,j)>=(m-1).and.pos(3,j)<m)hbbinnum(m)=hbbinnum(m)+1
                      enddo      
                      hbnum=hbnum+1
                    endif
                  else
                    v1(1)=pos(1,j)-pos(1,h2)
                    if(v1(1)>((hix-lox)/2))v1(1)=v1(1)-(hix-lox)
                    v1(2)=pos(2,j)-pos(2,h2)
                    if(v1(2)>((hiy-loy)/2))v1(2)=v1(2)-(hiy-loy)
                    v1(3)=pos(3,j)-pos(3,h2)
                    if(v1(3)>((hiz-loz)/2))v1(3)=v1(3)-(hiz-loz)
                    v2(1)=pos(1,k)-pos(1,h2)
                    if(v2(1)>((hix-lox)/2))v2(1)=v2(1)-(hix-lox)
                    v2(2)=pos(2,k)-pos(2,h2)
                    if(v2(2)>((hiy-loy)/2))v2(2)=v2(2)-(hiy-loy)
                    v2(3)=pos(3,k)-pos(3,h2)
                    if(v2(3)>((hiz-loz)/2))v2(3)=v2(3)-(hiz-loz)
                    v1m=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
                    v2m=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
                    v1v2=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
                    theta=acos( v1v2/(v1m*v2m) )
                    if(theta.ge.thetacut)then
                      do m=1,nbins
                        if(pos(3,j)>=(m-1).and.pos(3,j)<m)hbbinnum(m)=hbbinnum(m)+1
                      enddo
                      hbnum=hbnum+1
                    endif
                  endif
                endif
              endif
            enddo 
          endif
        endif
      enddo
    endif
  enddo
  !$OMP END PARALLEL DO
  print*,hbnum,h2onum
  write(*,'(F12.6,A12)'),dble(hbnum)/dble(h2onum),' HB/molecule'
  do k=1,nbins
    if(h2obinnum(k)>0)then
      print*,k,dble(hbbinnum(k))/dble(h2obinnum(k)),h2obinnum(k),hbbinnum(k)
    else
      print*,k,' 0.000000000000000E+000 ',h2obinnum(k),hbbinnum(k)
    endif
  enddo
  call cpu_time(finish)
  delta=secnds(t1)
  print '("Step ",I5," CPU Time = ",f12.6," seconds.")',i,finish-start
  print '("Step ",I5," Wall Time = ",f12.6," seconds.")',i,delta
enddo

end program hbonds


