double precision function q(vx,vy)
  implicit none                
  double precision :: vx,vy
  if(vy>=0) then
    q = acos(vx/sqrt((vx**2+vy**2)))
  else
    q= -acos(vx/sqrt((vx**2+vy**2)))
  end if
  return              
end function q 

double precision function EL(x,y,rx,ry,a,b,theta)
  implicit none                
  double precision :: x,y,rx,ry,a,b,theta
   EL = ((x-rx)*cos(theta)+(y-ry)*sin(theta))**2/a**2+(-(x-rx)*sin(theta)+(y-ry)*cos(theta))**2/b**2-1
  return
end function EL

double precision function US(x)
  implicit none                
  double precision :: x
   if(x>=0) then
     US = 1.0
   else 
     US = 0.0
   end if
  return
end function US


program main
  implicit none

  integer, parameter :: trial=10000, SP=16, step=10, initial=1000, add=1, cell=initial+add*int(trial/step)
  double precision :: q, EL, US
  integer :: access

  integer(8) i, j, n, k, Fi, slice, label(cell), nearest(cell), time0, time1, dtime
  integer(8)  near, countFa, countFr, j0(cell), nst, List(cell,cell)
  double precision DList(cell,cell), rate, xmax, ymax, r0x, r0y, time
  double precision pi, dt, dSP, theta0, q0, fl(cell), rxij, fad, frd, dv2, sqv
  double precision a1(cell), b1(cell), a2(cell), b2(cell)
  double precision EaPx, EaPy, ErPx1, ErPy1, dx, dy, Fax(cell), Fay(cell), Fa1x, Fa1y, Fr1x, Fr1y
  double precision Frx(cell), Fry(cell), ftheta(cell), fp
  double precision UStepSin, Friction, Friction2, Mass, maxa1, maxa2
  double precision gamma0, epsilon(cell), fa, fr, d(cell)
  double precision gamma1(cell), lambda, gamma2, w0, dphi,V2
  double precision rx(cell),ry(cell), vx(cell),vy(cell), phi(cell), psi(cell)
  double precision Lpair(3), rx2, ry2, dl(3), gamma1_Ctr, gamma1_KO, flattening_Ctr, flattening_KO
  double precision EP4(SP),Distance(cell,cell), SinPsi1, CosPsi1, SinEP(SP), CosEP(SP)
  character filename*128, path*128, ns*8
    
  integer, allocatable, dimension(:) :: seed
  integer :: nrand
  integer, parameter :: SEED_SELF = 1
  double precision :: Urandom1,Urandom2

  call random_seed(size=nrand)
  allocate(seed(nrand))
  seed = SEED_SELF
  call random_seed(put=seed)

! Start time
call system_clock(time0)

!Control : KO = rate : 1-rate
  rate = 0.0
  gamma1_Ctr = 0.005
  gamma1_KO  = 0.01
  flattening_Ctr = 0.75
  flattening_KO  = 0.75


!--------------------------------------------------- 
  pi = atan(1.0)*4.0
  dt = 1.0 
  dSP = dble(SP)
  slice = int(trial/2000)
  xmax=2500; ymax=200
  Lpair(1)=-1.0; Lpair(2)= 0.0; Lpair(3)= 1.0
!--------------------------------------------------- 
  a1 = 30.0
  b1 = 20.0
  gamma0 = 0.05
  epsilon = 0.15
  fa = 0.0075
  fr = 0.075
  d = 0.005
  lambda = 0.3
  gamma2 = 0.01
  w0 = 0.01

if(rate==1.0)then
    path = "Ctr_KO_1_0"
else if(rate==0.8)then
    path = "Ctr_KO_4_1"
else if(rate==0.5)then
    path = "Ctr_KO_1_1"
else if(rate==0.0)then
    path = "Ctr_KO_0_1"
end if

!label=0: KO, label=1: Control
label=0
do i=1, cell
    call random_number(Urandom1)
    if(1-rate > Urandom1)then
        label(i)=0
        gamma1(i)=gamma1_KO
		fl(i) = flattening_KO
    else
        label(i)=1
        gamma1(i)=gamma1_Ctr
		fl(i) = flattening_Ctr
    end if
end do

do i=1, cell
   a2(i) = 20.0
   b2(i) = a2(i)*(1.0-fl(i))
end do



rx=0.0; ry=0.0; vx=0.0; vy=0.0
do i=1,cell
    call random_number(Urandom1)
    call random_number(Urandom2)
        rx(i) = 20.0 + (xmax-40.0)*Urandom1
        ry(i) = 20.0 + (ymax-20.0)*Urandom2
end do

do i=1,cell
    call random_number(Urandom1)
      theta0=Urandom1*2*pi
    call random_number(Urandom1)
    call random_number(Urandom2)
      vx(i) = 0.1*Urandom1*Cos(theta0)
      vy(i) = 0.1*Urandom2*Sin(theta0)
end do

do k=1, SP
     SinEP(k) = sin((2*pi*dble(k-1))/dSP)
     CosEP(k) = cos((2*pi*dble(k-1))/dSP)
     EP4(k) = sin((4*pi*dble(k-1))/dSP)
end do

do i=1, cell
call random_number(Urandom1)
    psi(i)=Urandom1*pi*2
end do
phi=0.0


!---create a directory---
if(access( trim(path), " ")==2) then
  call system("mkdir " //trim(path))
end if

if(access( trim(path)//"\data", " ")==2) then
  call system("mkdir " //trim(path)//"\data")
end if

if(access( trim(path)//"\images", " ")==2) then
  call system("mkdir " //trim(path)//"\images")
end if


!-------------------------------------------------------------------------
 
  open (18, file=trim(path)//"/parameters.csv", status='replace') 
    write (18, *) 'trial',',', 'cell',',', 'step',',', 'slice'
    write (18, *) trial,',',cell,',',step,',',slice

    write (18,*)
	write (18,*) 'a1',',', 'b1',',', 'gamma0',',', 'fa',',', 'fr',',', 'd'
	write (18,*) a1(1),',',b1(1),',', gamma0,',',fa,',',fr,',',d(1)

	write (18,*)
	write (18,*) 'gamma1',',', 'fp',',', 'lambda',',', 'gamma2',',', 'w0'
	write (18,*)  gamma1(1),',',fp,',',lambda,',',gamma2,',',w0
  close (18)

!-------------------------------------------------------------------------

Fi=initial
    filename = trim(path)//"/data/00000000.txt"
    open (15, file=filename, status='replace') 
        do i=1, Fi
            write (15, *) rx(i),ry(i),psi(i),a2(i),b2(i), label(i)
        end do
    close (15)

Fax=0.0;Fay=0.0
Frx=0.0;Fry=0.0
ftheta=0.0

do n=1, trial

if(mod(n,step)==0) then
   Fi = Fi+add
end if

if(mod(n,slice)==0) then
    write(ns, '(i8.8)') int(n/slice)
    filename = trim(path)//"/data/"//ns//".txt"
    open (16, file=filename, status='replace') 
        do i=1, Fi
            write (16, *) rx(i),ry(i),psi(i),a2(i),b2(i),label(i)
        end do
    close (16)
end if


maxa1 = 2.0*maxval(a1)
maxa2 = 2.0*maxval(a2)

Distance=0.0
do i=1, Fi-1
    do j=i+1, Fi
        Distance(i,j) = Hypot(rx(j)-rx(i),ry(j)-ry(i))
        Distance(j,i) = Hypot(rx(i)-rx(j),ry(i)-ry(j))
    end do
end do

!---Interaction---
j0=0; List=0; DList=0.0
Fax=0.0; Fay=0.0; Frx=0.0; Fry=0.0

if(Fi>1) then
 do i=1, Fi-1
    do j=i+1, Fi
       rxij = rx(j)-rx(i)
	   dl=(/abs(rxij-xmax),abs(rxij),abs(rxij+xmax)/)
	   near = Minloc(dl,DIM=1)
	   rx2=rx(j)+Lpair(near)*xmax
	   ry2=ry(j)
	   countFa = 0
       countFr = 0
	   if(Distance(i,j)<=maxa1) then
          SinPsi1 = sin(psi(i))
          CosPsi1 = cos(psi(i))
          do k=1, SP
	         EaPx = rx(i) + a1(i)*CosPsi1*CosEP(k) - b1(i)*SinPsi1*SinEP(k)
             EaPy = ry(i) + a1(i)*SinPsi1*CosEP(k) + b1(i)*CosPsi1*SinEP(k)
             if( EL(EaPx,EaPy,rx2,ry2,a1(j),b1(j),psi(j))<0 ) then
	             countFa = 1
                 exit
             end if
          end do ! end of k-loop

		  if(countFa>0) then
		      j0(i) = j0(i) + 1
		      List(i,j0(i)) = j
		      j0(j) = j0(j) + 1
		      List(j,j0(j)) = i
			  DList(i, j0(i)) = Distance(i,j)
			  DList(j, j0(j)) = Distance(i,j)
              fad = fa/Distance(i,j)
	          Fa1x = fad*(rx(j)-rx(i))
		      Fa1y = fad*(ry(j)-ry(i))
              Fax(i) = Fax(i) + Fa1x
              Fay(i) = Fay(i) + Fa1y
              Fax(j) = Fax(j) - Fa1x
              Fay(j) = Fay(j) - Fa1y
		  end if
        end if

        if(Distance(i,j)<=maxa2) then
          SinPsi1 = sin(psi(i))
          CosPsi1 = cos(psi(i))
          do k=1, SP
              ErPx1 = rx(i) + a2(i)*CosPsi1*CosEP(k) - b2(i)*SinPsi1*SinEP(k)
              ErPy1 = ry(i) + a2(i)*SinPsi1*CosEP(k) + b2(i)*CosPsi1*SinEP(k)
              if( EL(ErPx1,ErPy1,rx2,ry2,a2(j),b2(j),psi(j))<0 ) then 
                  countFr = 1
                  exit
	          end if
          end do ! end of k-loop
          if(countFr>0) then
             frd = fr/Distance(i,j)
		     Fr1x = frd*(rx(i)-rx(j))
		     Fr1y = frd*(ry(i)-ry(j))
			 Frx(i) = Frx(i) + Fr1x
             Fry(i) = Fry(i) + Fr1y
   	         Frx(j) = Frx(j) - Fr1x
             Fry(j) = Fry(j) - Fr1y
		  end if

	  end if
  end do ! end of j-loop
end do   ! end of i-loop
end if

!---the nearest neighbor of cell-i---

do i=1, Fi
   if(j0(i)>0) then
      nst = Minloc(DList(i,1:j0(i)),DIM=1)
      nearest(i) = List(i,nst)
   else
      nearest(i) = 0
   end if
end do


!---time evolution of velocity, position and angle---
do i=1, Fi
    UStepSin = (1.0-US(sin(phi(i))+epsilon(i)))
    Friction = gamma0+(1.0-gamma0)*(UStepSin)
    Friction2 = 1.0-Friction
    Mass = 1.0-UStepSin

    sqv = sqrt(vx(i)**2+vy(i)**2)
    if(j0(i)>0 .and. sqv > 10E-20) then
        dv2 = d(i)/sqv
        dx = dv2*vx(i)
        dy = dv2*vy(i)
    else
        dx = 0.0
        dy = 0.0
    end if


    vx(i) = Friction2*vx(i) + Mass*(Fax(i) + Frx(i) + dx)
    vy(i) = Friction2*vy(i) + Mass*(Fay(i) + Fry(i) + dy)
    r0x=rx(i)
    r0y=ry(i)
    rx(i) = rx(i) + vx(i)
    ry(i) = ry(i) + vy(i)
    if(ry(i)<0.0)then
        ry(i) = r0y
    end if
    if(rx(i)>xmax)then
        rx(i) = rx(i) - xmax !mod(rx(i),xmax)
    elseif(rx(i)<0.0)then
        rx(i) = rx(i) + xmax
    end if

    if( sqrt(vx(i)**2+vy(i)**2) < 10E-16) then
       q0 = psi(i)
    else
       q0 = q(vx(i),vy(i))
    end if           
       psi(i) = psi(i) - gamma1(i)*sin(2*(psi(i)-q0))

    call random_number(Urandom1)
       if( Urandom1 < 0.5) then
          dphi = 0.0001
       else
          dphi = -0.0001
       end if

    V2=0.0
    if(j0(i)>0) then
      if(label(i)==1 .and. label(nearest(i))==1)then
	      V2=(lambda*sin(phi(i)-phi(nearest(i))) + w0)
      else
          V2=0.0
      end if
    else
      V2=0.0
    end if
      phi(i) = phi(i) - gamma2*sin(phi(i)) + V2 + dphi

end do !i-loop
end do !n-loop

! Finish time
      call system_clock(time1, dtime)
! Output time
      time = 1d0*(time1-time0)/dtime
      print'(" Time (s)",f12.4)', time

  stop
end program main
