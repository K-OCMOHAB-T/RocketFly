program jopa
    
    double precision Mol, g0, Rk, R
    parameter ( Mol = 0.02896442d+0, g0 = 9.8066d+0, Rk = 8.31432d+0, R = 6356767.0d+0 )
    
    double precision B, gp, gu, gd
    parameter ( B = Mol*g0/Rk, gp = 0.00005324d+0, gu = 20.046796, gd = Mol/Rk )
    
    double precision P(11), Hd(10)
    
    double precision Sr, h0
    parameter ( Sr = 26.364d+0, h0 = 91.0d+0 )
    
    double precision t0, t01, t1, tb, t2, t23, t3
    parameter ( t0 = 21.5d+0, t01 = 61.0d+0, t1 = 117.8d+0, tb = 157.48d+0, t2 = 287.3d+0, t23 = 396.0d+0, t3d = 524.96d+0, t3 = 528.26d+0 )
    !parameter ( t0 = 21.5d+0, t01 = 61.0d+0, t1 = 117.8d+0, tb = 157.48d+0, t2 = 287.3d+0, t23 = 396.0d+0, t3 = 528.26d+0 )
    
    double precision mf1, me1, mf2, me2, mf3, me3, mp, msas, mb
    parameter ( mf1 = 43400.0d+0, me1 = 3810.0d+0, mf2 = 99500.0d+0, me2 = 6550.0d+0, mf3 = 25300.0d+0, me3 = 2410.0d+0, mp = 7300.0d+0, msas = 2300.0d+0, mb = 2200.0d+0 )
    !parameter ( mf1 = 43400.0d+0, me1 = 3810.0d+0, mf2 = 99500.0d+0, me2 = 6550.0d+0, mf3 = 25300.0d+0, me3 = 2410.0d+0, mp = 7300.0d+0, msas = 2300.0d+0, mb = 2200.0d+0 )
    
    
    double precision fs1, fv1, fs2, fv2, fv3
    parameter ( fs1 = 838500.0d+0, fv1 = 1021300.0d+0, fs2 = 792480.0d+0, fv2 = 990180.0d+0, fv3 = 297930.0d+0 )
    
    double precision gradtorad, radtograd
    parameter ( gradtorad = 0.0174532925d+0, radtograd = 57.2957795131d+0 )
    double precision v01, v1, v2, v23, v3
    parameter ( v01 = gradtorad*60.3d+0, v1 = gradtorad*32.7d+0, v2 = gradtorad*19.5d+0, v23 = gradtorad*11.9d+0, v3 = gradtorad*1.0d+0 )
    !parameter ( v01 = gradtorad*60.5d+0, v1 = gradtorad*32.5d+0, v2 = gradtorad*19.6d+0, v23 = gradtorad*12.0d+0, v3 = gradtorad*1.0d+0 )
    
    type :: PDUtype
        double precision P
        double precision D
        double precision Usound
        logical Doff
    end type
    
    double precision time, dt, dt2, dt6, ux, uy, x, y, h
    double precision k_ux, k_uy, k_x, k_y, kux, kuy, kx, ky, kxmem, kymem
    double precision dux, duy, dx, dy
    integer i
    type(PDUtype) pdu
    
    Hd(1) = 11000.0d+0
    Hd(2) = 20000.0d+0
    Hd(3) = 32000.0d+0
    Hd(4) = 47000.0d+0
    Hd(5) = 51000.0d+0
    Hd(6) = 71000.0d+0
    Hd(7) = 85000.0d+0
    Hd(8) = 94000.0d+0
    Hd(9) = 102450.0d+0
    Hd(10) = 117777.0d+0
    
    P(1) = 101325.0d+0
    P(11) = 0.00266618d+0
    

    dt = 0.01d+0
    
    dt2 = dt/2.0d+0
    dt6 = dt/6.0d+0
    time = 0.0d+0
    
    do i = 1, 9
        pdu = PDUcalc(time, Hd(i))
        P(i+1) = pdu%P
        write(*,*) i, Hd(i), P(i+1)
    enddo
    
    
    
    x = 0.0d+0
    y = h0
    ux = 0.0d+0
    uy = 0.0d+0
    
    open(11, file="loghjhjj.txt", status="replace")
    write(11,'(a15,1X,a15,1X,a15,1X,a15)') 'time', 'F/mass', 'R/mass', 'g'
    
    open(12, file="log.txt", status="replace")
    
    write(12,'(a15,1X,a15,1X,a15,1X,a15,1X,a15,1X,a15,1X,a15,1X,a15)') 'time', 'x', 'y', 'l', 'h', 'ux', 'uy', 'A/g'
    write(12,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') time, x, y, lreport(x,y), hreport(x,y), ux, uy
    i = 1
    do while (time .le. t3+1.0)
    
        ! 1 step
        kx = ux
        ky = uy
        h = x*x
        h = sqrt(h + (y+R)*(y+R)) - R
        pdu = PDUcalc(time, h)
        kux = Ax(time, x, y, ux, uy, h, pdu)
        kuy = Ay(time, x, y, ux, uy, h, pdu)
    
        dx = kx
        dy = ky
        dux = kux
        duy = kuy
    
        k_x = kx
        k_y = ky
        k_ux = kux
        k_uy = kuy
        
        ! 2 step
        kx = ux + dt2*k_ux
        ky = uy + dt2*k_uy
        kxmem = x + dt2*k_x
        kymem = y + dt2*k_y
        h = kxmem*kxmem
        h = sqrt(h + (kymem+R)*(kymem+R)) - R
        pdu = PDUcalc(time + dt2, h)
        kux = Ax(time + dt2, kxmem, kymem, kx, ky, h, pdu)
        kuy = Ay(time + dt2, kxmem, kymem, kx, ky, h, pdu)
    
        dx = dx + kx*2.0d+0
        dy = dy + ky*2.0d+0
        dux = dux + kux*2.0d+0
        duy = duy + kuy*2.0d+0
    
        k_x = kx
        k_y = ky
        k_ux = kux
        k_uy = kuy
        
        ! 3 step
        kx = ux + dt2*k_ux
        ky = uy + dt2*k_uy
        kxmem = x + dt2*k_x
        kymem = y + dt2*k_y
        h = kxmem*kxmem
        h = sqrt(h + (kymem+R)*(kymem+R)) - R
        pdu = PDUcalc(time + dt2, h)
        kux = Ax(time + dt2, kxmem, kymem, kx, ky, h, pdu)
        kuy = Ay(time + dt2, kxmem, kymem, kx, ky, h, pdu)
    
        dx = dx + kx*2.0d+0
        dy = dy + ky*2.0d+0
        dux = dux + kux*2.0d+0
        duy = duy + kuy*2.0d+0
    
        k_x = kx
        k_y = ky
        k_ux = kux
        k_uy = kuy
        
        ! 4 step
        kx = ux + dt*k_ux
        ky = uy + dt*k_uy
        kxmem = x + dt*k_x
        kymem = y + dt*k_y
        h = kxmem*kxmem
        h = sqrt(h + (kymem+R)*(kymem+R)) - R
        pdu = PDUcalc(time + dt, h)
        kux = Ax(time + dt, kxmem, kymem, kx, ky, h, pdu)
        kuy = Ay(time + dt, kxmem, kymem, kx, ky, h, pdu)
    
        dx = dx + kx
        dy = dy + ky
        dux = dux + kux
        duy = duy + kuy
    
        x = x + dt6*dx
        y = y + dt6*dy
        ux = ux + dt6*dux
        uy = uy + dt6*duy
        
        time = time + dt
        
        if ((mod(i, 100) .eq. 0) .or. (time .ge. t3)) then
            call Areport(time, x, y, ux, uy)
            write(12,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') time, x, y, lreport(x,y), hreport(x,y), ux, uy, sqrt(kux*kux + kuy*kuy)/g0, angle(time)*radtograd
        endif
        i = i + 1
    enddo
    
    close(12)
    
    write(*,*) cos(1.570796330d+0)
    
    

    call system('pause')
   
   
    
    contains
    
    double precision function hreport(x, y)
    implicit none
    double precision x, y  
        hreport = x*x
        hreport = sqrt(hreport + (y+R)*(y+R)) - R
    end function
    
    double precision function lreport(x, y)
    implicit none
    double precision x, y
        lreport = R*atan(x/(y+R))
    end function
    

    double precision function Ax(t, x, y, ux, uy, h, pdu)
    implicit none
    double precision t, x, y, ux, uy, h
    type(PDUtype) pdu
    double precision curvate
    
        curvate = atan(x/(y+R))
    
        Ax = (Fforce(t, pdu%P) - Rforce(t, ux, uy, pdu))/mass(t) * cos(angle(t) - curvate) - grav(h)*sin(curvate)
        
    end function
    
    double precision function Ay(t, x, y, ux, uy, h, pdu)
    implicit none
    double precision t, x, y, ux, uy, h
    type(PDUtype) pdu
    double precision curvate
    
        curvate = atan(x/(y+R))
    
        Ay = (Fforce(t, pdu%P) - Rforce(t, ux, uy, pdu))/mass(t) * sin(angle(t) - curvate) - grav(h)*cos(curvate)
        
        !write(11,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') t, Fforce(t, pdu%P)/mass(t), Rforce(t, ux, uy, pdu)/mass(t), grav(h)
        
    end function
    
    subroutine Areport(t, x, y, ux, uy)
    implicit none
    double precision t, x, y, ux, uy
    type(PDUtype) pdu
    double precision curvate, h
        h = x*x
        h = sqrt(h + (y+R)*(y+R)) - R
        
        pdu = PDUcalc(t, h)
    
        curvate = atan(x/(y+R))
        
        write(11,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') t, Fforce(t, pdu%P), Rforce(t, ux, uy, pdu), grav(h)*mass(t)
        
    end subroutine
    
    
    double precision function Rforce(t, ux, uy, pdu)
    implicit none
    double precision t, ux, uy
    type(PDUtype) pdu
    
    double precision density, usound, cx_m, u, u2
    
        u2 = ux*ux + uy*uy
        u = sqrt(u2)
    
        if (t .le. t1) then
            density = pdu%D
            usound = pdu%Usound
            cx_m = Cx(u/usound)
            
            Rforce = cx_m*density*u2*Sr*0.5d+0
            
            !write(11,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') t, density, usound, cx_m, Rforce
        else
            cx_m = Cx_d_s_2(u)
            
            Rforce = cx_m*u2
        endif
    end function
    
    double precision function Fforce(t, press)
    implicit none
    double precision t, press
    
        if (t .le. t1) then
            Fforce = fv2 - press/P(1) * (fv2 - fs2)
            Fforce = Fforce + (fv1 - press/P(1) * (fv1 - fs1))*4.0d+0
        else if (t .le. t2) then
            Fforce = fv2 - press/P(1) * (fv2 - fs2)
        else if (t .le. t3d) then
            Fforce = fv3
        else
            Fforce = 0.0d+0
        endif
    end function
    
    
    double precision function mass(t)
    implicit none
    double precision t
    
        if (t .le. t1) then
            mass = mb + msas + mp
            mass = mass + mf3
            mass = mass + (me2*t - mf2*(t-t2))/t2
            mass = mass + 4.0d+0*(me1*t - mf1*(t-t1))/t1
        else if (t .le. tb) then
            mass = mb + mp
            mass = mass + mf3
            mass = mass + (me2*t - mf2*(t-t2))/t2
        else if (t .le. t2) then
            mass = mp
            mass = mass + mf3
            mass = mass + (me2*t - mf2*(t-t2))/t2
        else if (t .le. t3d) then
            mass = mp + (me3*(t-t2) - mf3*(t-t3d))/(t3d-t2)
        else if (t .le. t3) then
            mass = mp + me3
        else
            mass = mp
        endif
    end function
    
    double precision function grav(h)
    implicit none
    double precision h
        grav = h/R + 1.0d+0
        grav = g0/grav/grav
    end function
    
    double precision function angle(t)
    implicit none
    double precision t
    double precision pi2, mem
    pi2 = 1.570796330d+0
    
        if (t .le. t0) then
            angle = pi2
        else if (t .le. t01) then
            mem = (t01 - t)/(t01 - t0)
            angle = v01 + (pi2 - v01)*mem*mem
        else if (t .le. t1) then
            mem = (t1 - t)/(t1 - t01)
            angle = v1 + (v01 - v1)*mem*mem
        !else if (t .le. t1) then
        !    angle = v01 + (v1 - v01)*(t - t01)/(t1 - t01)
        !else if (t .le. t1) then
        !    mem = (t1 - t)/(t1 - t0)
        !    angle = v1 + (pi2 - v1)*mem*mem
        !else if (t .le. t2) then
        !    angle = v1 + (v2 - v1)*(t - t1)/(t2 - t1)
        else if (t .le. t2) then
            mem = (t2 - t)/(t2 - t1)
            angle = v2 + (v1 - v2)*mem*mem
        else if (t .le. t23) then
            angle = v2 + (v23 - v2)*(t - t2)/(t23 - t2)
        else if (t .le. t3) then
            angle = v23 + (v3 - v23)*(t - t23)/(t3 - t23)
        else
            angle = 0.0d+0
        endif
    end function
    
    
    double precision function Cx(max)
    implicit none
    double precision max
    
    double precision maxpol(6)
    integer i
    
        maxpol(1) = max
        do i = 2, 6
            maxpol(i) = maxpol(i-1)*max
        enddo
        
        if (max .le. 0.6d+0) then
            Cx = -26.3381239d+0*maxpol(6)
            Cx = Cx + 46.1635160d+0*maxpol(5)
            Cx = Cx - 28.9618844d+0*maxpol(4)
            Cx = Cx + 7.4424147d+0*maxpol(3)
            Cx = Cx - 0.4767275d+0*maxpol(2)
            Cx = Cx - 0.1051217d+0*maxpol(1)
            Cx = Cx + 0.3901999
        else if (max .le. 1.15d+0) then
            Cx = -33.0301991d+0*maxpol(5)
            Cx = Cx + 118.6021743d+0*maxpol(4)
            Cx = Cx - 162.0887050d+0*maxpol(3)
            Cx = Cx + 105.3394111d+0*maxpol(2)
            Cx = Cx - 32.2073362d+0*maxpol(1)
            Cx = Cx + 3.9811664d+0
        else if (max .le. 1.579311d+0) then
            Cx = -426.617373d+0
            Cx = Cx + 1550.244014d+0*maxpol(1)
            Cx = Cx - 2233.689317d+0*maxpol(2)
            Cx = Cx + 1598.580663d+0*maxpol(3)
            Cx = Cx - 568.510363d+0*maxpol(4)
            Cx = Cx + 80.401721d+0*maxpol(5)
        else if (max .le. 6.0d+0) then
            Cx = 1.4835155878d+0
            Cx = Cx - 0.7722261483d+0*maxpol(1)
            Cx = Cx + 0.0934728027d+0*maxpol(2)
            Cx = Cx + 0.0579681350d+0*maxpol(3)
            Cx = Cx - 0.0245109115d+0*maxpol(4)
            Cx = Cx + 0.0036081641d+0*maxpol(5)
            Cx = Cx - 0.0001896281d+0*maxpol(6)
        else
            Cx = 0.18d+0
        endif
    
    end function
    
    double precision function Cx_d_s_2(u)
    implicit none
    double precision u
        Cx_d_s_2 = 64.17641628d+0*exp(-0.00497522d+0*u)
    end function
   
    type(PDUtype) function PDUcalc(t,h)
    implicit none
    double precision h, t
    double precision alfa, temp1, temp
    
        if (h .le. Hd(1)) then
            alfa = -0.0065d+0
            temp1 = 288.15d+0
            temp = alfa*h + temp1
            
            PDUcalc%P = Pressure( h, 0.0d+0, alfa, temp1, temp, 1, .false. )

        else if (h .le. Hd(2)) then
            alfa = 0.0d+0
            temp1 = 216.65d+0
            temp = temp1
            
            PDUcalc%P = Pressure( h, Hd(1), alfa, temp1, temp, 2, .true. )

        else if (h .le. Hd(3)) then
            alfa = 0.001d+0
            temp1 = 216.65d+0
            temp = alfa*(h-Hd(2)) + temp1
            
            PDUcalc%P = Pressure( h, Hd(2), alfa, temp1, temp, 3, .false. )

        else if (h .le. Hd(4)) then
            alfa = 0.0028d+0
            temp1 = 228.65d+0
            temp = alfa*(h-Hd(3)) + temp1
            
            PDUcalc%P = Pressure( h, Hd(3), alfa, temp1, temp, 4, .false. )

        else if (h .le. Hd(5)) then
            alfa = 0.0d+0
            temp1 = 270.65d+0
            temp = temp1
            
            PDUcalc%P = Pressure( h, Hd(4), alfa, temp1, temp, 5, .true. )

        else if (h .le. Hd(6)) then
            alfa = -0.0028d+0
            temp1 = 270.65d+0
            temp = alfa*(h-Hd(5)) + temp1
            
            PDUcalc%P = Pressure( h, Hd(5), alfa, temp1, temp, 6, .false. )

        else if (h .le. Hd(7)) then
            alfa = -0.002d+0
            temp1 = 214.65d+0
            temp = alfa*(h-Hd(6)) + temp1
            
            PDUcalc%P = Pressure( h, Hd(6), alfa, temp1, temp, 7, .false. )

        else if (h .le. Hd(8)) then
            alfa = 0.0d+0
            temp1 = 186.65d+0
            temp = alfa*(h-Hd(7)) + temp1
            
            PDUcalc%P = Pressure( h, Hd(7), alfa, temp1, temp, 8, .true. )

        else if (h .le. Hd(9)) then
            alfa = 0.002d+0
            temp1 = 186.65d+0
            temp = alfa*(h-Hd(8)) + temp1
            
            PDUcalc%P = Pressure( h, Hd(8), alfa, temp1, temp, 9, .false. )

        else
            alfa = 0.0085d+0
            temp1 = 203.55d+0
            temp = temp1
            
            PDUcalc%P = Pressure( h, Hd(9), alfa, temp1, temp, 10, .false. )

        endif

        
        if (t .lt. t1) then
            PDUcalc%D = gd*PDUcalc%P/temp
            PDUcalc%Usound = gu*sqrt(temp)
            
            !write(11,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') t, PDUcalc%P, temp, PDUcalc%D
        endif
        
    end function
    
    double precision function Pressure( h, h1, a, temp1, temp, i, simple )
    implicit none
    double precision h, h1, a, temp1, temp
    integer i
    logical simple
    
    double precision mem1, mem2, mem3
    
    !write(*,*) h, Hd(10), P(i)
        if (h .lt. Hd(10)) then
            mem1 = h/R + 1.0d+0
            mem2 = h1/R + 1.0d+0
            
            if (simple) then
                Pressure = P(i)*exp(-B*(h-h1)/temp1/mem1/mem2)
            else
                mem3 = a*(R+h1)
                
                Pressure = P(i)*exp(-B*(h-h1)/(temp1-mem3)/mem1/mem2)
                
                mem3 = a*R
                mem3 = mem2 - temp1/mem3
                mem3 = a*mem3*mem3
                
                Pressure = Pressure*exp(B/mem3*log(temp1/temp*mem1/mem2))
            endif
        else
            Pressure = P(11)*exp(-gp*(h-Hd(10)))
        endif
    end function
    
    
end