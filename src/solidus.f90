!==========================================================
! Solidus Subroutine
!==========================================================
! Four different solidus may be chosed from for melt production models
subroutine solidus
use arrays
use params
implicit none
integer :: i,j
double precision :: c1,c2,c3,c4,Tolm
double precision :: stressone,xPs,dPdTm,tint,Tsm,deltam
double precision, external :: stressI

!$OMP Parallel private(i,j,stressone,xPs,dPdTm,tint,Tsm,deltam)

!==========================================================
! Calculate liquidus
!==========================================================

! McKenzie and Bickle liquidus
!

if (isolidus.eq.4) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tliq(j,i) =  1736.2 + 4.343*stressone + 180.*atan(stressone/2.2169)
        end do
    end do
    !$OMP end do

! Katz 2003 DEFAULT
else

    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tliq(j,i) =  1780. + 45.*stressone - 2.*(stressone**2)
        end do
    end do 
    !$OMP end do
end if

if (icrust_melt.eq.1) then
! Liquidus for upper and lower crust
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tliquc(j,i) =  (1262. + 0.09*stressone*1000.) - 273.
            tliqlc(j,i) =  (1423 +  0.105*stressone*1000.) - 273.
        end do
    end do 
    !$OMP end do
endif

!==========================================================
! Calculate Solidus
!==========================================================

! Upper and lower crust solidus (Gerya p372)
 if (icrust_melt.eq.1) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            ! Change pressure to GPa
            stressone = -stressI(j,i)/1.e6
            if (stressone.le.1.2) then
                tsoluc(j,i) = (889.+17900./(stressone+54.)+20200./(stressone+54.)**2.)-273.
            else
                tsoluc(j,i) = (831.+0.06*stressone)-273.
            end if
            if (stressone.le.1.6) then
                tsollc(j,i) = (973.+70400./(stressone+354.)+77800000./(stressone+354.)**2.)-273. 
            else
                tsollc(j,i) = (935.+0.0035*stressone+0.0000062*stressone*stressone)-273.
            end if
        end do
    end do
    !$OMP end do
endif

! Katz et al (2003) DEFAULT
if (isolidus.eq.1) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tsol(j,i) = (-5.104)*(stressone**2) + 132.899*(stressone) + 1085.7
        end do
    end do
    !$OMP end do

! Hirschmann et al (2000) & Hirschmann et al (2009) 
else if (isolidus.eq.2) then
    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            ! Change pressure to GPa
            stressone = -stressI(j,i)/1.e9
            if (stressone.le.10.) then
                tsol(j,i) = -5.140*stressone*stressone + 132.899*stressone + 1120.661 
            elseif (stressone.le.23.5) then
                tsol(j,i) = -1.092*((stressone-10.)**2) + 32.390*(stressone-10.) + 1935.
            else
                tsol(j,i) = 26.35*(stressone-23.5) + 2175
            end if
        end do
    end do
!$OMP end do
! Herzberg et al (2000) 
else if (isolidus.eq.3) then

!$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            if (stressone.le.2.7) then
                tsol(j,i) = 132*(stressone) + 1102
            else
                tsol(j,i) = 390*(log(stressone)) - 5.7*(stressone) + 1086
            end if
        end do
    end do
!$OMP end do

! McKenzie & Bickle (1988) Solidus 
else if (isolidus.eq.4) then
    c1 = 1100.
    c2 = 136.
    c3 = 0.0004968
    c4 = 0.012
    Tolm = 0.01

    !$OMP do
    do j = 1,nz-1
        do i = 1,nx-1
            stressone = -stressI(j,i)/1.e9
            tsol(j,i) = 1000. 
            xPs = ((tsol(j,i)-c1)/c2)+(c3*exp(c4*(tsol(j,i)-c1)))
            dPdTm = (c3*c4*exp(c4*(tsol(j,i)-c1)) + (1/c2))
            tint = xPs - (dPdTm*tsol(j,i))
            Tsm = (stressone-tint)/dPdTm
            deltam = abs(Tsm-tsol(j,i))
            tsol(j,i) = Tsm
            
            do while (deltam.ge.Tolm)
                xPs = ((tsol(j,i)-c1)/c2)+(c3*exp(c4*(tsol(j,i)-c1)))
                dPdTm= (c3*c4*exp(c4*(tsol(j,i)-c1)))+(1/c2)
                tint = xPs - (dPdTm*tsol(j,i))
                Tsm = (stressone-tint)/dPdTm
                deltam = abs(Tsm-tsol(j,i))
                tsol(j,i) = Tsm
            end do
        end do
    end do
    !$OMP end do
endif
!$OMP end parallel 

return
end subroutine solidus
