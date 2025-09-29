!===========================================
! Melt Production Code adapted from Bown & White (1995) and McKenzie &
! Bickle (1988)
!===========================================

subroutine melt_prod
use arrays
use params
use phases
implicit none
integer :: i, j, m, iblk, jblk
double precision :: xpress, dumFF, dumP, DTsol, dumsol, dumF, dumT, arbar
double precision :: xpTpP, hP, xpTpF, AvDF, AvDT
double precision :: xDF(4), xDT(4)
double precision, external :: CalcF, CalcT, DeltaTwat, calcsolid, stressI
double precision, dimension(4) :: dumC

!dumC_local(:) = dumC(:)  ! Workaround to avoid shared array reference in OMP clause
dumC(1) = 0.0 !Runge Kutta
dumC(2) = 0.5 !Runge Kutta
dumC(3) = 0.5 !Runge Kutta
dumC(4) = 1.0 !Runge Kutta

!===========================================
! Nodes to Elements, T and P calc, and Initialization
!===========================================
!$OMP parallel do collapse(2) schedule(static,1) default(none) &
!$OMP private(i,j,xpress,dumFF,dumP,DTsol,dumsol,dumF,dumT,arbar,xpTpP,xpTpF,hP,AvDF,AvDT) &
!$OMP shared(temp, avT, phase_ratio, Eff_melt,Melt_extract, fmelt, nloop, nx, nz, watercont, xmodalcpx)
do j = 1,nz-1-10
    do i = 1+5,nx-1-5
        avT(j,i) = 0.25*(temp(j,i) + temp(j,i+1) + temp(j+1,i) + temp(j+1,i+1))
        !  Lithostatic pressure
        xpress = -stressI(j,i)/1.e9 ! Pressure in GPa   
        
        !===========================
        ! Initialize
        !===========================
        if(nloop.eq.1) then
            Eff_melt(j,i) = 0.
            Melt_extract(j,i) = 0.
            if (phase_ratio(kmant2,j,i).gt.0.95) then  
                ! Calculate initial melt
                Eff_melt(j,i) = CalcF(xpress,avT(j,i),watercont,xmodalcpx)
                fmelt(j,i) = Eff_melt(j,i) ! Only for the first time step. 
                Melt_extract(j,i) = Eff_melt(j,i)
            end if
      end if
    end do
end do
!$OMP end parallel do


!=================================
! Calculate changes in melt and temps
!=================================
do iblk = 1,2
do jblk = 1,2
!$OMP parallel do collapse(2) schedule(static,1) default(none) &
!$OMP private(i,j,m,xpress,dumFF,dumP,DTsol,dumsol,hP,dumF,dumT,arbar,xpTpP,xpTpF,AvDF,AvDT,xDF,xDT) &
!$OMP shared(jblk, iblk, nloop, phase_ratio, fmagma, tsol, zpressm, Eff_melt,Melt_extract, fmelt, avT, temp, &
!$OMP deltaTLH, watercont, xmodalcpx,dumC, nx,nz)
do j = jblk,nz-35,2
    do i = iblk+10,nx-11,2
        xpress = -stressI(j,i)/1.e9 ! Pressure in GPa   
    
        if (phase_ratio(kmant2,j,i).gt.0.95.and.nloop.ne.1.and.fmagma(j,i).eq.0.) then         
            dumFF = 0.
            dumP = xpress
            DTsol = deltaTwat(watercont,dumP,dumFF)
            dumsol = tsol(j,i) 

            if (avT(j,i).gt.(dumsol-DTsol)) then
                hP = (zpressm(j,i)/1.e9) - xpress !Amount of decompression in GPa

                if (hP.gt.1e-12) then
                    do m = 1,4 !4th order solver
                        dumP = (zpressm(j,i)/1.e9) - hP*dumC(m)

                        if (m == 1) then
                            dumF = Melt_extract(j,i)
                            dumT = avT(j,i)
                        else
                            dumF = Melt_extract(j,i) - dumC(m)*hP*xDF(m-1)
                            dumT = avT(j,i) - dumC(m)*hP*xDT(m-1)
                        end if

                        arbar = dumF*(68./2.9) + (1.-dumF)*(40./3.3)
                        ! Eqn 22 Katz 2003
                        xpTpP = (calcT(dumP+1.e-14,watercont,xmodalcpx,dumF) &
                                -calcT(dumP-1.e-14,watercont,xmodalcpx,dumF)) &
                                / 2.e-14

                        dumFF = 0.
                        DTsol = deltaTwat(watercont,dumP,dumFF)
                        dumsol = calcsolid(dumP)

                        if (dumT.gt.(dumsol-DTsol)) then
                            dumF = max(dumF,1.e-9) !Work out from a reasonable value
                            ! Eqn 21 Katz 2003
                            xpTpF = (calcT(dumP,watercont,xmodalcpx,dumF+1.e-14) &
                                -calcT(dumP,watercont,xmodalcpx,dumF-1.e-14)) &
                                / 2.e-14
                            ! Melt Derivative - dFdP - Eqn 20 Katz2003
                                xDF(m) = (-1000.0*xpTpP/(dumT+273.0) + arbar)/(1000.0*xpTpF/ &
                                (dumT+273.0)+ 300.0)
                        else
                            xDF(m) = 0.
                        end if

                        xDT(m) = (dumT+273.0)*(arbar - 300.0*xDF(m))/1000.0 
                    end do

                    AvDF = (xDF(1)/6. + xDF(2)/3. + xDF(3)/3. + xDF(4)/6.)
                    AvDT = (xDT(1)/6. + xDT(2)/3. + xDT(3)/3. + xDT(4)/6.)

                    if (AvDF.le.0.) then
                        !$OMP ATOMIC UPDATE
                        Melt_extract(j,i) = Melt_extract(j,i) + abs(AvDF * hP) 
                        Eff_melt(j,i) = Eff_melt(j,i) + abs(AvDF*hP) ! Keep 5% melt in asthenosphere
                        !$OMP ATOMIC WRITE
                        fmelt(j,i) = abs(AvDF*hP) ! 95% is migrating?

                        if (Eff_melt(j,i).ge.1.) then
                            Eff_melt(j,i) = 0.99
                            fmelt(j,i) = 0.
                            AvDT = 0.0
                        end if

                        ! Proportionally distribute the element change in temperature to the nodes
                        if (AvDF.ne.0.) then
                            deltaTLH(j,i)= AvDT*hP
                            avT(j,i) = avT(j,i) - deltaTLH(j,i)

                            temp(j,i)     = temp(j,i)     - deltaTLH(j,i) * (temp(j,i)     / avT(j,i) / 4.)
                            temp(j,i+1)   = temp(j,i+1)   - deltaTLH(j,i) * (temp(j,i+1)   / avT(j,i) / 4.)
                            temp(j+1,i)   = temp(j+1,i)   - deltaTLH(j,i) * (temp(j+1,i)   / avT(j,i) / 4.)
                            temp(j+1,i+1) = temp(j+1,i+1) - deltaTLH(j,i) * (temp(j+1,i+1) / avT(j,i) / 4.)
                        end if
                    end if
                end if
            end if
        end if
    end do
end do
!$OMP end parallel do
end do
end do

return
end subroutine melt_prod
