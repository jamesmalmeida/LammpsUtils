program perfildensidade
implicit none
integer                                     :: ii,jj,kk,ll,ccount,nnumber,ntype,hcount,ocount,sicount
integer                                     :: natoms,nbins,dummyint,nconf,counter
double precision                            :: dr,raio,low,high,lo,hi,d,r,numerador,denominador
double precision                            :: Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,x,y,z
double precision                            :: xh,xl,yh,yl,dl,dh,lo1,hi1,zl,zh,lox,loy,hix,hiy
double precision,dimension(:),allocatable   :: xx,yy,zz
integer,dimension(:),allocatable            :: nntype,nnnumber
character(50)                               :: dummychar,inputfile,outputfile
INTEGER                                     :: atomsinside
integer                                     :: dummy
double precision,parameter                  :: pi=3.1415926536

integer                                     :: nmol,nmoltype,nspecies
integer,dimension(:),allocatable            :: nnmol,moltype(:)
integer,allocatable                         :: num(:,:),nr(:,:,:)
double precision,allocatable                :: atmmass(:),massnr(:,:),dens(:,:,:),densm(:)


!-----Read Input File-----!
!        open(9,file='in.dens-2D',status='old')
!        read(9,*)inputfile
!        read(9,*)outputfile
!        read(9,*)nbins
!        read(9,*)raio
!        read(9,*)nconf
!        read(9,*)jump
!        close(9)
!-------------------------!
nbins=218
nconf=800
nspecies=38
nmoltype=11

open(10,file='positions.lammpstrj',status='old')

read(10,*)dummychar
read(10,*)dummyint
read(10,*)dummychar
read(10,*)natoms
read(10,*)dummychar
read(10,*)lo,hi
read(10,*)lo,hi
read(10,*)lo1,hi1
read(10,*)dummychar

do ii=1,natoms
 read(10,*)nnumber,ntype,dummy,x,y,z!,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz
enddo

allocate(xx(1:natoms),yy(1:natoms),zz(1:natoms),nntype(1:natoms),nnnumber(1:natoms))
allocate(nnmol(1:natoms),num(1:nmoltype,1:nspecies),nr(1:nmoltype,1:nspecies,1:nbins),atmmass(1:nspecies),moltype(1:natoms))
allocate(massnr(1:nmoltype,1:nbins),dens(1:nconf,1:nmoltype,nbins),densm(1:nmoltype))

dr=(hi1-lo1)/nbins !raio/nbins
dh=0.0
dl=100.0

atmmass(1)=(1E24/6.02E23)*12.011
atmmass(2)=(1E24/6.02E23)*15.9994
atmmass(3)=(1E24/6.02E23)*14
atmmass(4)=(1E24/6.02E23)*12.011
atmmass(5)=(1E24/6.02E23)*1.0079
atmmass(6)=(1E24/6.02E23)*12.011
atmmass(7)=(1E24/6.02E23)*1.0079
atmmass(8)=(1E24/6.02E23)*12.011
atmmass(9)=(1E24/6.02E23)*1.0079
atmmass(10)=(1E24/6.02E23)*12.011
atmmass(11)=(1E24/6.02E23)*1.0079
atmmass(12)=(1E24/6.02E23)*12.011
atmmass(13)=(1E24/6.02E23)*1.0079
atmmass(14)=(1E24/6.02E23)*15.9994
atmmass(15)=(1E24/6.02E23)*1.0079
atmmass(16)=(1E24/6.02E23)*35.5
atmmass(17)=(1E24/6.02E23)*22.98
atmmass(18)=(1E24/6.02E23)*32
atmmass(19)=(1E24/6.02E23)*15.9994
atmmass(20)=(1E24/6.02E23)*24.3
atmmass(21)=(1E24/6.02E23)*40.078
atmmass(22)=(1E24/6.02E23)*39.1
atmmass(23)=(1E24/6.02E23)*12.011
atmmass(24)=(1E24/6.02E23)*32
atmmass(25)=(1E24/6.02E23)*15.9994
atmmass(26)=(1E24/6.02E23)*1.0079
atmmass(27)=(1E24/6.02E23)*15.9994
atmmass(28)=(1E24/6.02E23)*1.0079
atmmass(29)=(1E24/6.02E23)*15.9994
atmmass(30)=(1E24/6.02E23)*15.9994
atmmass(31)=(1E24/6.02E23)*12.011
atmmass(32)=(1E24/6.02E23)*12.011
atmmass(33)=(1E24/6.02E23)*15.9994
atmmass(34)=(1E24/6.02E23)*15.9994
atmmass(35)=(1E24/6.02E23)*1.0079
atmmass(36)=(1E24/6.02E23)*32
atmmass(37)=(1E24/6.02E23)*12.011
atmmass(38)=(1E24/6.02E23)*1.0079

rewind(10)

do ll=1,nconf
  read(10,*)dummychar
  read(10,*)dummyint
  read(10,*)dummychar
  read(10,*)natoms
  read(10,*)dummychar
  read(10,*)lox,hix
  read(10,*)loy,hiy
  read(10,*)lo,hi
  read(10,*)dummychar

  do ii=1,natoms
    read(10,*)nnumber,ntype,nmol,x,y,z!,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz
    !if(x>=hix/2)x=x-(hix-lox)
    !if(y>=hiy/2)y=y-(hiy-loy)
    xx(ii)=x !-((hi-lo)/2)!-((xh-xl)/2)
    yy(ii)=y !-((hi-lo)/2)!-((yh-yl)/2)
    zz(ii)=z !-((hi-lo)/2)
    nnnumber(ii)=nnumber
    nntype(ii)=ntype
    nnmol(ii)=nmol
    if(ntype==14)moltype(nmol)=1  !H2O
    if(ntype==34)moltype(nmol)=2  !H3O
    if(ntype==31)moltype(nmol)=3  !Prot
    if(ntype==32)moltype(nmol)=4  !DeProt
    if(ntype==16)moltype(nmol)=5  !Cl
    if(ntype==17)moltype(nmol)=6  !Na
    if(ntype==18)moltype(nmol)=7  !S
    if(ntype==19)moltype(nmol)=8  !O
    if(ntype==20)moltype(nmol)=9  !Mg
    if(ntype==8) moltype(nmol)=10 !C
    if(ntype==22)moltype(nmol)=11 !K
  enddo

  !!!!DEBUG1!!!!
  !  do ii=1,natoms
  !    write(*,*)nnnumber(ii),nntype(ii),nnmol(ii),moltype(nnmol(ii))
  !  enddo
  !  stop
  !!!!DEBUG1-!!!! MOLTYPES ARE BEING ATRIBUTED CORRECTLY

  atomsinside=0

  do jj=1,nbins
    num=0
    low=lo+(jj-1)*dr
    high=lo+jj*dr

    do ii=1,natoms
      !----- Boundary Conditions -----!
        !xx(ii)=xx(ii)-nint(xx(ii)/hi)*hi
        !yy(ii)=yy(ii)-nint(yy(ii)/hi)*hi
      !----- Boundary Conditions -----!
      d=zz(ii)
      if ((d.ge.low).and.(d.le.high))then
        do kk=1,nmoltype
          if(moltype(nnmol(ii))==kk)then
            select case (nntype(ii))
             case(1)
               num(kk,1)=num(kk,1)+1
             case(2)
               num(kk,2)=num(kk,2)+1
             case(3)
               num(kk,3)=num(kk,3)+1
             case(4)
               num(kk,4)=num(kk,4)+1
             case(5)
               num(kk,5)=num(kk,5)+1
             case(6)
               num(kk,6)=num(kk,6)+1
             case(7)
               num(kk,7)=num(kk,7)+1
             case(8)
               num(kk,8)=num(kk,8)+1
             case(9)
               num(kk,9)=num(kk,9)+1
             case(10)
               num(kk,10)=num(kk,10)+1
             case(11)
               num(kk,11)=num(kk,11)+1
             case(12)
               num(kk,12)=num(kk,12)+1
             case(13)
               num(kk,13)=num(kk,13)+1
             case(14)
               num(kk,14)=num(kk,14)+1
             case(15)
               num(kk,15)=num(kk,15)+1
             case(16)
               num(kk,16)=num(kk,16)+1
             case(17)
               num(kk,17)=num(kk,17)+1
             case(18)
               num(kk,18)=num(kk,18)+1
             case(19)
               num(kk,19)=num(kk,19)+1
             case(20)
               num(kk,20)=num(kk,20)+1
             case(21)
               num(kk,21)=num(kk,21)+1
             case(22)
               num(kk,22)=num(kk,22)+1
             case(23)
               num(kk,23)=num(kk,23)+1
             case(24)
               num(kk,24)=num(kk,24)+1
             case(25)
               num(kk,25)=num(kk,25)+1
             case(26)
               num(kk,26)=num(kk,26)+1
             case(27)
               num(kk,27)=num(kk,27)+1
             case(28)
               num(kk,28)=num(kk,28)+1
             case(29)
               num(kk,29)=num(kk,29)+1
             case(30)
               num(kk,30)=num(kk,30)+1
             case(31)
               num(kk,31)=num(kk,31)+1
             case(32)
               num(kk,32)=num(kk,32)+1
             case(33)
               num(kk,33)=num(kk,33)+1
             case(34)
               num(kk,34)=num(kk,34)+1
             case(35)
               num(kk,35)=num(kk,35)+1
             case(36)
               num(kk,36)=num(kk,36)+1
             case(37)
               num(kk,37)=num(kk,37)+1
             case(38)
               num(kk,38)=num(kk,38)+1
            end select
          endif
        enddo        
      endif
    enddo
    do kk=1,nmoltype
      do ii=1,nspecies 
        nr(kk,ii,jj)=num(kk,ii)
      enddo
    enddo
  enddo

  !!!!DEBUG2!!!!
  !do jj=1,nbins
  !  !do kk=1,nmoltype
  !    do ii=1,nspecies
  !      !write(*,*)kk,ii,jj,nr(kk,ii,jj)
  !      write(*,*)ii,jj,nr(:,ii,jj)
  !    enddo
  !  !enddo
  !enddo
  !stop
  !!!!DEBUG2-!!!! Seems to be counting OK, different number for each molecule type

  massnr=0

  do jj=1,nbins
    do kk=1,nmoltype
      do ii=1,nspecies
        massnr(kk,jj)=massnr(kk,jj)+nr(kk,ii,jj)*atmmass(ii)
      enddo
    enddo

    denominador=dr*(hix-lox)*(hiy-loy)

    do kk=1,nmoltype
      dens(ll,kk,jj)=massnr(kk,jj)/denominador 
    enddo
  enddo
  write(*,*)ll
enddo

!!!!DEBUG3!!!!
!do jj=1,nbin
!    write(*,'(I7,4F12.4)')jj,massnr(:,jj)
!enddo
!stop
!!!!DEBUG3-!!!! Seems perfect

!!!!DEBUG4!!!!
!!do ll=1,nconf
!!  !do kk=1,nmoltype
!!    do jj=1,nbins
!!      !write(*,*)jj,kk,ii,dens(ll,kk,jj)
!!      write(*,'(2I7,4F15.4)')ll,jj,dens(ll,:,jj)
!!    enddo
!!  !enddo
!!enddo
!!stop
!!!!DEBUG4-!!!! Seems perfect 

close(10)
open(13,file='dens-2D.dat')

do jj=1,nbins
  densm=0
  r=(jj-1)*dr

  do ll=1,nconf
    do kk=1,nmoltype
      densm(kk)=densm(kk)+dens(ll,kk,jj)
    enddo
  enddo

  do kk=1,nmoltype
    densm(kk)=densm(kk)/(nconf)
  enddo
  write(13,"(20F12.3)")r,densm(:)
  !write(*,"(17F12.3)")r,densm(:)
enddo

close(13)
close(14)
!deallocate(xx(1:natoms),yy(1:natoms),zz(1:natoms),nntype(1:natoms),nnnumber(1:natoms))
!deallocate(nnmol(1:natoms),num(1:nmoltype,1:nspecies),nr(1:nmoltype,1:nspecies,1:nbins),atmmass(1:nspecies),moltype(1:natoms))
!deallocate(massnr(1:nspecies,1:nbins),dens(1:nconf,1:nmoltype,nbins),densm(1:nmoltype))

end program perfildensidade

