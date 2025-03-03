subroutine marker2elem
  !$ACC routine gang
  use myrandom_mod
  use marker_data
  use arrays
  use params
  use phases
  implicit none

  integer :: kph(1), i, j, k, kinc, inc, iseed, icount
  double precision :: x1, x2, y1, y2, xx, yy, r1, r2
  double precision :: kn, xdepth, testu
  integer :: mohod
  double precision :: critmeltfrac(nz-1,nx-1), tmpr, critmelt, testa
  integer :: ibug,imagcolumn(nx-1),jmagcolumn(nx-1),jmagtop,magtopcount,imag_min,imag_max, knn
  
  !character*200 msg

  ! Interpolate marker properties into elements
  ! Find the triangle in which each marker belongs

  !$OMP parallel do private(kinc,r1,r2,x1,y1,x2,y2,xx,yy,inc,icount)
  !$ACC loop collapse(2) gang vector
  do i = 1 , nx-1
      do j = 1 , nz-1
          iseed = nloop + i + j
          kinc = nmark_elem(j,i)
          icount = 0

          !  if there are too few markers in the element, create a new one
          !  with age 0 (similar to initial marker)
          !if(kinc.le.4) then
          !    write(msg,*) 'marker2elem: , create a new marker in the element (i,j))', i, j
          !    call SysMsg(msg)
          !endif

          do while (kinc.le.4)
              call myrandom(iseed, r1)
              call myrandom(iseed, r2)

              ! (x1, y1) and (x2, y2)
              x1 = cord(j  ,i,1)*(1-r1) + cord(j  ,i+1,1)*r1
              y1 = cord(j  ,i,2)*(1-r1) + cord(j  ,i+1,2)*r1
              x2 = cord(j+1,i,1)*(1-r1) + cord(j+1,i+1,1)*r1
              y2 = cord(j+1,i,2)*(1-r1) + cord(j+1,i+1,2)*r1

              ! connect
              ! (this point is not uniformly distributed within the element area
              ! and is biased against the thicker side of the element, but this
              ! point is almost gauranteed to be inside the element)
              xx = x1*(1-r2) + x2*r2
              yy = y1*(1-r2) + y2*r2

              call add_marker(xx, yy, iphase(j,i), zpressm(j,i), Eff_melt(j,i), 0.d0, j, i, inc)
              icount = icount + 1
              if(icount > 100) stop 133
              if(inc.le.0) cycle

              kinc = kinc + 1
              zpresscounter(j,i) = zpresscounter(j,i) + zpressm(j,i)
              Emeltcounter(j,i) = Emeltcounter(j,i) + Eff_melt(j,i)
          enddo

          Eff_melt(j,i) = Emeltcounter(j,i)/float(kinc)
          zpressm(j,i) = zpresscounter(j,i)/float(kinc)

          call count_phase_ratio(j,i)

      enddo
  enddo
  !$OMP end parallel do

  ! Find the Moho

jmoho(:) = 31 
do i = 1, nx-1
    kn = 0
    do j = 1, nz-1
        xdepth = cord(j+1,i,2)
        testu = sum(phase_ratio(mantle_phases,j,i))
        if ( testu > 0.5 .and. kn == 0) then
            if (phase_ratio(1,j,i) > 0.1.or.fmelt(j,i)>0.) then
                mohod = j - 1
            else
                mohod = j
            endif
            kn = kn + 1
        endif
    enddo
    jmoho(i) = mohod
    if (jmoho(i).le.8) jmoho(i) = 8
enddo

imagtop = int((nx-1)/2)
imagcolumn(:) = 0
jmagcolumn(:) = nz-1
ibug = 0

critmeltfrac = 0.

! find top magma chambers
do i = 1, nx-1
    k = jmoho(i)    
    do j = k , nz-1
        tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
        if (tmpr.ge.1100.) then
            critmeltfrac(j,i) = Eff_melt(j,i)
        endif
    enddo
enddo

critmelt = maxval(critmeltfrac)
do i = 1, nx-1
    k = jmoho(i)
    do j = k , nz-1
        if (critmeltfrac(j,i) .ge. critmelt) then
            jmagcolumn(i) = j
            imagcolumn(i) = i
            ibug = ibug + 1
            cycle
        endif
    enddo
enddo

jmagtop = minval(jmagcolumn)
print *, imagtop,jmagtop,ibug
magtopcount = 0
! find top magma chambers
if (ibug .gt. 0 ) then
    do i = 1, nx-1
        if (jmagcolumn(i) .eq. jmagtop.and.magtopcount.eq.0) then
            magtopcount = magtopcount + 1
            imag_min = imagcolumn(i)
        endif
        if (jmagcolumn(i) .eq. jmagtop.and. magtopcount .gt. 0) then
            magtopcount = magtopcount + 1
            imag_max = imagcolumn(i)
        endif
    enddo
    if (imag_max.eq.0) then
        imagtop = imag_min
    else
        imagtop = int(0.5*(imag_max+imag_min))
    endif
    print *, imagtop
endif

! The bottom of the magma need to be 15 km above the max melt fraction ot We produce too much melt
!  Find the basement below the extrusives and sediments
ibasement = 2! first guess below the extrusives
knn = 0
do j = 1, nz-1
    testa = sum(phase_ratio(surface_phases,j,imagtop))
    if ( testa> 0.5d0 .and. knn == 0) then
        ibasement = j 
        knn = knn + 1
    endif
enddo

if (ibasement.le.2) then
    ibasement = 2
else
    ibasement = max(2,ibasement)
endif

!   !$OMP parallel do
!   !$ACC loop auto
!   do i = 1, nx-1
!       jmoho(i) = nz-1
!       do j = 1, nz-1
!           if (sum(phase_ratio(mantle_phases,j,i)) > 0.5d0) then
!               jmoho(i) = j
!               exit
!           endif
!       enddo
!       !print *, i, jmoho(i)
!   enddo
!   !$OMP end parallel do

  return
end subroutine marker2elem
