program jopa
    
    double precision Mol, g0, Rk, R, pi2
    parameter ( Mol = 0.02896442d+0, g0 = 9.8066d+0, Rk = 8.31432d+0, R = 6356767.0d+0, pi2 = 1.570796330d+0 )
    
    double precision B, gp, gu, gd
    parameter ( B = Mol*g0/Rk, gp = 0.00005324d+0, gu = 20.0393403730729d+0, gd = Mol/Rk )
    
    double precision P(11), Hd(10)
    
    double precision Sr1, Sr2, h0
    parameter ( Sr1 = 26.364d+0, Sr2 = 7.072d+0, h0 = 91.0d+0 )
    
    double precision t0, t01, t1, tb, t2, t23, t3
    parameter ( t0 = 21.5d+0, t01 = 61.0d+0, t1 = 117.8d+0, tb = 157.48d+0, t2 = 287.3d+0, t23 = 396.0d+0, t3d = 524.96d+0, t3 = 528.26d+0 )
    
    double precision mf1, me1, mf2, me2, mf3, me3, mp, msas, mb
    parameter ( mf1 = 43400.0d+0, me1 = 3810.0d+0, mf2 = 99500.0d+0, me2 = 6550.0d+0, mf3 = 25300.0d+0, me3 = 2410.0d+0, mp = 7300.0d+0, msas = 1950.0d+0, mb = 2200.0d+0 )
    
    
    double precision fs1, fv1, fs2, fv2, fv3
    parameter ( fs1 = 838500.0d+0, fv1 = 1021300.0d+0, fs2 = 792480.0d+0, fv2 = 990180.0d+0, fv3 = 297930.0d+0 )
    !parameter ( fs1 = 838500.0d+0, fv1 = 1021300.0d+0, fs2 = 792480.0d+0, fv2 = 990180.0d+0, fv3 = 297930.0d+0 )
    
    double precision gradtorad, radtograd
    parameter ( gradtorad = 0.0174532925d+0, radtograd = 57.2957795131d+0 )
    double precision v01, v1, v2, v23, v3
    parameter ( v01 = gradtorad*60.3d+0, v1 = gradtorad*32.7d+0, v2 = gradtorad*19.5d+0, v23 = gradtorad*11.9d+0, v3 = gradtorad*1.0d+0 )
    
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
    
    double precision maxx, alfaa
    
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
    write(11,'(a15,1X,a15,1X,a15,1X,a15,1X,a15,1X,a15)') 'time', 'alfa', 'F/mass', 'Rx/mass', 'Ry/mass', 'g'
    
    
    !maxx = 0.0d+0
    !alfaa = 7.5d+0*gradtorad
    !do while (maxx .le. 5.3d+0)
    !    write(11,'(f15.5,1X,f15.5,1X,f15.5)') maxx, Cxs1(maxx, alfaa), Cys1(maxx, alfaa)
    !    maxx = maxx + 0.01d+0
    !enddo
    !do while (maxx .le. 16.0d+0)
    !    write(11,'(f15.5,1X,f15.5,1X,f15.5)') maxx, Cxs2(maxx, alfaa), Cys2(maxx, alfaa)
    !    maxx = maxx + 0.01d+0
    !enddo
    !call system('pause')
    
    
    
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
    double precision curvate, cosp, sinp, alfa, psi
    
        curvate = atan(x/(y+R))
        
        sinp = sqrt(ux*ux + uy*uy)
        if (sinp .ne. 0.0d+0) then 
            cosp = ux/sinp
            sinp = uy/sinp
        else
            cosp = 0.0d+0
            sinp = 0.0d+0
        endif
        if (cosp .gt. 1.0d+0) cosp = 1.0d+0
        if (sinp .gt. 1.0d+0) sinp = 1.0d+0
        
        if (ux .ne. 0.0d+0) then 
            psi = atan(uy/ux)
        else
            psi = pi2
        endif
        
        alfa = angle(t) - curvate - psi

        Ax = (Fforce(t, pdu%P) * cos(angle(t) - curvate) - Rxforce(t, ux, uy, pdu, alfa) * cosp - Ryforce(t, ux, uy, pdu, alfa) * sinp)/mass(t) - grav(h)*sin(curvate)
        
    end function
    
    double precision function Ay(t, x, y, ux, uy, h, pdu)
    implicit none
    double precision t, x, y, ux, uy, h
    type(PDUtype) pdu
    double precision curvate, cosp, sinp, alfa, psi
    
        curvate = atan(x/(y+R))
        
        sinp = sqrt(ux*ux + uy*uy)
        if (sinp .ne. 0.0d+0) then 
            cosp = ux/sinp
            sinp = uy/sinp
        else
            cosp = 0.0d+0
            sinp = 0.0d+0
        endif
        if (cosp .gt. 1.0d+0) cosp = 1.0d+0
        if (sinp .gt. 1.0d+0) sinp = 1.0d+0
        
        if (ux .ne. 0.0d+0) then 
            psi = atan(uy/ux)
        else
            psi = pi2
        endif
        
        alfa = angle(t) - curvate - psi
    
        Ay = (Fforce(t, pdu%P) * sin(angle(t) - curvate) - Rxforce(t, ux, uy, pdu, alfa) * sinp + Ryforce(t, ux, uy, pdu, alfa) * cosp)/mass(t) - grav(h)*cos(curvate)
        
    end function
    
    subroutine Areport(t, x, y, ux, uy)
    implicit none
    double precision t, x, y, ux, uy
    type(PDUtype) pdu
    double precision curvate, h, psi
        h = x*x
        h = sqrt(h + (y+R)*(y+R)) - R
        
        pdu = PDUcalc(t, h)
    
        curvate = atan(x/(y+R))
        
        if (ux .ne. 0.0d+0) then 
            psi = atan(uy/ux)
        else
            psi = pi2
        endif
        
        write(11,'(f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5,1X,f15.5)') t, (angle(t) - curvate - psi)*radtograd,Fforce(t, pdu%P), Rxforce(t, ux, uy, pdu, angle(t) - curvate - psi), Ryforce(t, ux, uy, pdu, angle(t) - curvate - psi), grav(h)
        
    end subroutine
    
    
    double precision function Rxforce(t, ux, uy, pdu, alfa)
    implicit none
    double precision t, ux, uy, alfa
    type(PDUtype) pdu
    
    double precision density, usound, cx, u, u2
    
        u2 = ux*ux + uy*uy
        u = sqrt(u2)
    
        density = pdu%D
        usound = pdu%Usound
        
        if (t .le. t1) then
            cx = Cxs1(u/usound, alfa)
            Rxforce = cx*density*u2*Sr1*0.5d+0
        else
            cx = Cxs2(u/usound, alfa)
            Rxforce = cx*density*u2*Sr2*0.5d+0
        endif
    end function
    
    double precision function Ryforce(t, ux, uy, pdu, alfa)
    implicit none
    double precision t, ux, uy, alfa
    type(PDUtype) pdu
    
    double precision density, usound, cy, u, u2
    
        u2 = ux*ux + uy*uy
        u = sqrt(u2)
    
        density = pdu%D
        usound = pdu%Usound
        
        if (t .le. t1) then
            cy = Cys1(u/usound, alfa)
            Ryforce = cy*density*u2*Sr1*0.5d+0
        else
            cy = Cys2(u/usound, alfa)
            Ryforce = cy*density*u2*Sr2*0.5d+0
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
    double precision mem
    
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
    
    double precision function Cys1(max, alfa)
    implicit none
    double precision max, alfa
    
    double precision al, acof, bcof
    double precision maxpol(6)
    integer i
    
        maxpol(1) = max
        do i = 2, 6
            maxpol(i) = maxpol(i-1)*max
        enddo
        
        ! Cy = a * alfa*alfa + b * alfa
        ! a
        if (max .le. 0.2d+0) then
            acof = 0.20847579d+0*max
            acof = acof - 1.86527273d+0
        else if (max .le. 0.8d+0) then
            acof = -2.31501290d+0*maxpol(2)
            acof = acof + 1.82898482d+0*maxpol(1)
            acof = acof - 2.09677402d+0
        else if (max .le. 1.15d+0) then
            acof = -30.07464381d+0*maxpol(2)
            acof = acof + 60.89962573d+0*maxpol(1)
            acof = acof - 31.58073696d+0
        else if (max .le. 5.290801d+0) then
            acof = 17.58882875d+0
            acof = acof - 51.31951242d+0*maxpol(1)
            acof = acof + 49.70608299d+0*maxpol(2)
            acof = acof - 21.77501810d+0*maxpol(3)
            acof = acof + 4.93703377d+0*maxpol(4)
            acof = acof - 0.55947624d+0*maxpol(5)
            acof = acof + 0.02483795d+0*maxpol(6)
        else
            acof = 11.44268707d+0
            acof = acof - 0.94515280d+0*max
        endif
        
        ! b
        if (max .le. 0.9d+0) then
            bcof = 1.37006770d+0*maxpol(2)
            bcof = bcof - 0.51129728d+0*maxpol(1)
            bcof = bcof + 3.54991076d+0
        else if (max .le. 1.15d+0) then
            bcof = -103.33882948d+0
            bcof = bcof + 321.71850827d+0*maxpol(1)
            bcof = bcof - 317.87893179d+0*maxpol(2)
            bcof = bcof + 103.55054523d+0*maxpol(3)
        else if (max .le. 1.3d+0) then
            bcof = 4.84648732d+0
            bcof = bcof - 0.97085573d+0*max
        else if (max .le. 2.05301d+0) then
            bcof = 2.46199276d+0
            bcof = bcof + 1.65293266d+0*maxpol(1)
            bcof = bcof - 0.60735524d+0*maxpol(2)
        else if (max .le. 4.44305d+0) then
            bcof = 10.47172958d+0
            bcof = bcof - 9.13831047d+0*maxpol(1)
            bcof = bcof + 4.29395325d+0*maxpol(2)
            bcof = bcof - 0.89066441d+0*maxpol(3)
            bcof = bcof + 0.06718546d+0*maxpol(4)
        else
            bcof = 2.83478397d+0
            bcof = bcof - 0.03078069d+0*max
        endif
        
        !write(11,'(f15.5,1X,f15.5,1X,f15.5)') max, acof, bcof
        
        if (alfa .lt. 0.0d+0) then
            al = -alfa
            Cys1 = -acof*al*al - bcof*al
        else
            al = alfa
            Cys1 = acof*al*al + bcof*al
        endif
        
    
    end function
    
    double precision function Cys2(max, alfa)
    implicit none
    double precision max, alfa
    
    double precision dcx, al, al2, acof, bcof
        
        if (alfa .lt. 0.0d+0) then
            al = -alfa
        else
            al = alfa
        endif
        al2 = al*al
        
        if (max .le. 6.248604d+0) then
            acof = -2.72930222d+0*al2
            acof = acof + 0.20436333d+0*al
            bcof = 32.45447226d+0*al2
            bcof = bcof - 0.48460404d+0*al
        
            Cys2 = bcof + acof*max
        else if (max .le. 15.0d+0) then
            acof = -1.09901134d+0*al2
            acof = acof + 0.01354598d+0*al
            bcof = 22.52382088d+0*al2
            bcof = bcof + 0.67772891d+0*al
        
            Cys2 = bcof + acof*max
        else
            acof = -1.09901134d+0*al2
            acof = acof + 0.01354598d+0*al
            bcof = 22.52382088d+0*al2
            bcof = bcof + 0.67772891d+0*al
            
            Cys2 = bcof + acof*15.0d+0
        endif
        
        if (alfa .lt. 0.0d+0) Cys2 = -Cys2
    
    end function
    
    double precision function Cxs1(max, alfa)
    implicit none
    double precision max, alfa
    
    double precision cx0, al, acof, bcof
    double precision maxpol(6)
    integer i
    
        maxpol(1) = max
        do i = 2, 6
            maxpol(i) = maxpol(i-1)*max
        enddo
        
        ! Cx при нулевом угле атаки
        cx0 = Cx0s1(max, maxpol)
        
        if (alfa .lt. 0.0d+0) then
            al = -alfa
        else
            al = alfa
        endif
        
        ! dCx = a * alfa*alfa + b * alfa
        ! a
        if (max .le. 0.2d+0) then
            acof = 0.40485636d+0*max
            acof = acof + 3.32811318d+0
        else if (max .le. 0.8d+0) then
            acof = -2.00060513d+0*maxpol(2)
            acof = acof + 1.80527995d+0*maxpol(1)
            acof = acof + 3.12805266d+0
        else if (max .le. 1.0d+0) then
            acof = 86.68253534d+0*maxpol(2)
            acof = acof - 153.75776231d+0*maxpol(1)
            acof = acof + 70.82127657d+0
        else if (max .le. 2.05301d+0) then
            acof = 42.27362164d+0
            acof = acof - 79.47929612d+0*maxpol(1)
            acof = acof + 51.48700699d+0*maxpol(2)
            acof = acof - 10.54951235d+0*maxpol(3)
        else if (max .le. 5.290801d+0) then
            acof = 46.18230094d+0
            acof = acof - 49.43708713d+0*maxpol(1)
            acof = acof + 21.10417725d+0*maxpol(2)
            acof = acof - 3.85886874d+0*maxpol(3)
            acof = acof + 0.25776934d+0*maxpol(4)
        else
            acof = 3.23749409d+0
            acof = acof + 0.47026298d+0*max
        endif
        
        ! b
        if (max .le. 0.2d+0) then
            bcof = 0.17217694d+0*max
            bcof = bcof - 0.00583712d+0
        else if (max .le. 0.9d+0) then
            bcof = 3.51925670d+0*maxpol(3)
            bcof = bcof - 4.23247376d+0*maxpol(2)
            bcof = bcof + 1.76239845d+0*maxpol(1)
            bcof = bcof - 0.18273652d+0
        else if (max .le. 1.3d+0) then
            bcof = 0.75097877d+0
            bcof = bcof - 0.24125531d+0*max
        else if (max .le. 2.669855d+0) then
            bcof = 3.54496320d+0
            bcof = bcof - 3.91431231d+0*maxpol(1)
            bcof = bcof + 1.35467786d+0*maxpol(2)
            bcof = bcof - 0.13702744d+0*maxpol(3)
        else if (max .le. 5.290801d+0) then
            bcof = 3.65034204d+0
            bcof = bcof - 2.73342622d+0*maxpol(1)
            bcof = bcof + 0.68225454d+0*maxpol(2)
            bcof = bcof - 0.05630134d+0*maxpol(3)
        else
            bcof = 0.27195714d+0
            bcof = bcof - 0.05911962d+0*max
        endif
        
        Cxs1 = cx0 + acof*al*al + bcof*al
    
    end function
    
    double precision function Cxs2(max, alfa)
    implicit none
    double precision max, alfa
    
    double precision cx0, dcx, al, al2, acof, bcof, ccof
    double precision max2
    
        max2 = max*max
        
        ! Cx при нулевом угле атаки
        cx0 = Cx0s2(max)
        
        if (alfa .lt. 0.0d+0) then
            al = -alfa
        else
            al = alfa
        endif
        al2 = al*al
        
        ! dCx
        if (max .le. 5.329115d+0) then
            ccof = 5.09177268d+0*al2
            ccof = ccof + 0.30887437d+0*al
            
            dcx = ccof
        else if (max .le. 6.248604d+0) then  
            acof = -4.78874255d+0*al2
            acof = acof + 1.05711951d+0*al
            bcof = 57.06854655d+0*al2
            bcof = bcof - 12.34169099d+0*al
            ccof = -163.03535115d+0*al2
            ccof = ccof + 36.05753707d+0*al
            
            dcx = ccof + bcof*max + acof*max2
        else if (max .le. 15.0d+0) then 
            acof = 1.47106164d+0*al2
            acof = acof - 0.06849292d+0*al
            bcof = -2.36174018d+0*al2
            bcof = bcof + 0.58868354d+0*al
            
            dcx = bcof + acof*max
        else
            acof = 1.47106164d+0*al2
            acof = acof - 0.06849292d+0*al
            bcof = -2.36174018d+0*al2
            bcof = bcof + 0.58868354d+0*al
            
            dcx = bcof + acof*15.0d+0
        endif
        
        Cxs2 = cx0 + dcx
    
    end function
    
    
    double precision function Cx0s1(max, maxpol)
    implicit none
    double precision max, maxpol(*)
        
        if (max .le. 0.6d+0) then
            Cx0s1 = -26.3381239d+0*maxpol(6)
            Cx0s1 = Cx0s1 + 46.1635160d+0*maxpol(5)
            Cx0s1 = Cx0s1 - 28.9618844d+0*maxpol(4)
            Cx0s1 = Cx0s1 + 7.4424147d+0*maxpol(3)
            Cx0s1 = Cx0s1 - 0.4767275d+0*maxpol(2)
            Cx0s1 = Cx0s1 - 0.1051217d+0*maxpol(1)
            Cx0s1 = Cx0s1 + 0.3901999
        else if (max .le. 1.15d+0) then
            Cx0s1 = -33.0301991d+0*maxpol(5)
            Cx0s1 = Cx0s1 + 118.6021743d+0*maxpol(4)
            Cx0s1 = Cx0s1 - 162.0887050d+0*maxpol(3)
            Cx0s1 = Cx0s1 + 105.3394111d+0*maxpol(2)
            Cx0s1 = Cx0s1 - 32.2073362d+0*maxpol(1)
            Cx0s1 = Cx0s1 + 3.9811664d+0
        else if (max .le. 1.579311d+0) then
            Cx0s1 = -426.617373d+0
            Cx0s1 = Cx0s1 + 1550.244014d+0*maxpol(1)
            Cx0s1 = Cx0s1 - 2233.689317d+0*maxpol(2)
            Cx0s1 = Cx0s1 + 1598.580663d+0*maxpol(3)
            Cx0s1 = Cx0s1 - 568.510363d+0*maxpol(4)
            Cx0s1 = Cx0s1 + 80.401721d+0*maxpol(5)
        else if (max .le. 6.0d+0) then
            Cx0s1 = 1.4835155878d+0
            Cx0s1 = Cx0s1 - 0.7722261483d+0*maxpol(1)
            Cx0s1 = Cx0s1 + 0.0934728027d+0*maxpol(2)
            Cx0s1 = Cx0s1 + 0.0579681350d+0*maxpol(3)
            Cx0s1 = Cx0s1 - 0.0245109115d+0*maxpol(4)
            Cx0s1 = Cx0s1 + 0.0036081641d+0*maxpol(5)
            Cx0s1 = Cx0s1 - 0.0001896281d+0*maxpol(6)
        else
            Cx0s1 = 0.18d+0
        endif
    
    end function
    
    double precision function Cx0s2(max)
    implicit none
    double precision max
        if (max .le. 6.248604d+0) then
            Cx0s2 = 0.57059573d+0
            Cx0s2 = Cx0s2 - 0.01890772d+0*max
            Cx0s2 = Cx0s2 + 0.00824439d+0*max*max
        else if (max .le. 15.0d+0) then
            Cx0s2 = 0.29397706d+0
            Cx0s2 = Cx0s2 + 0.07680994d+0*max
        else
            Cx0s2 = 1.446126d+0
        endif
    end function
    
    !double precision function Cx_d_s_2(u)
    !implicit none
    !double precision u
    !    Cx_d_s_2 = 64.17641628d+0*exp(-0.00497522d+0*u)
    !end function
   
    type(PDUtype) function PDUcalc(t,h)
    implicit none
    double precision h, t
    double precision alfa, temp1, temp, usound, density
    
        ! Temperature + Pressure
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
        
        ! Usound
        if (h .lt. Hd(8)) then
            usound = gu*sqrt(temp)
        else
            usound = -0.00016326d+0 * h
            usound = usound + 289.251d+0
        endif
            

        ! Density
        if (t .lt. t1) then
            density = gd*PDUcalc%P/temp
        else
            if (h .lt. 65000.0d+0) then
                density = 0.56476857d+0 * exp(-0.00012667d+0 * h)
            else
                density = 1.84584296d+0 * exp(-0.00014498d+0 * h)
            endif
        endif
           
        PDUcalc%D = density
        PDUcalc%Usound = usound
        
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