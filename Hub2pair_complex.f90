subroutine Hub2Pair(x,y,h,v,w,NAO, &
    H00,H11,H11,H20,H40,H31,H22,Hb22,H13,H04)

Use Precision
Implicit None
Integer                        :: NAO
Complex(Kind=pr), Intent(in)   :: x(NAO), y(NAO)               ! x(NAO) and y(NAO) are the arrays for u and v respectively
Complex(Kind=pr), Intent(in)   :: h(NAO), v(NAO,NAO), w(NAO,NAO) ! h(NAO) is the one electron integrals, v,w are the two electron integrals
Complex(Kind=pr), Intent(out)  :: H00
Complex(Kind=pr), Intent(out)  :: H02(NAO), H11(NAO), H20(NAO)
Complex(Kind=pr), Intent(out)  :: H04(NAO,NAO), H13(NAO,NAO), H22(NAO,NAO)
Complex(Kind=pr), Intent(out)  :: H40(NAO,NAO), H31(NAO,NAO), Hb22(NAO,NAO)  ! Hb22 is the \tilde{H}^{22}

! To get the equations refer to my overleaf document

! Initialize H11 to zero
do p=1, NAO
    H11(p) = 0.0
end do

! First sum over p (single sums)
do p = 1, NAO
    H11(p) = H11(p) - ( &
        - conjg(x(p)) * h(p) * x(p) + (conjg(y(p)))**2 * v(p,p) * (y(p))**2 + conjg(y(p)) * h(p) * y(p) &
    )
end do

! Double sum over p,q
do p = 1, NAO
    do q = 1, NAO
        H11(p) = H11(p) - ( &
            conjg(x(p)) * conjg(y(p)) * v(q,p) * x(q) * y(q) &
            - conjg(x(p)) * conjg(y(q)) * w(p,q) * x(p) * y(q) &
            + conjg(x(q)) * conjg(y(q)) * v(p,q) * x(p) * y(p) &
            + conjg(y(p)) * conjg(y(q)) * w(p,q) * y(p) * y(q) &
        )
    end do
end do


do p=1,NAO
    do q=1,NAO
        Hb22(p,q) = 0.0  ! initialize complex 0.0  
    end do
end do  

do p = 1, NAO
    do q = 1, NAO
        Hb22(p,q) = Hb22(p,q) + &
            (conjg(x(q)))**2 * v(p,q) * x(p)**2 + &
            2*conjg(x(q)) * conjg(y(p)) * w(p,q) * x(p) * y(q) + &
            (conjg(y(p)))**2 * v(q,p) * y(q)**2
    end do
end do


do p=1,NAO
    do q=1,NAO
        H13(p,q) = 0.0  ! initialize complex 0.0  
    end do
end do 


 
do p = 1, NAO
    do q = 1, NAO
        H13(p,q) = H13(p,q) + &
            conjg(x(q)) * conjg(y(p)) * w(p,q) * x(p) * x(q) - &
            conjg(x(q)) * conjg(y(q)) * v(p,q) * x(p)**2 + &
            (conjg(y(p)))**2 * v(q,p) * x(q) * y(q) - &
            conjg(y(p)) * conjg(y(q)) * w(p,q) * x(p) * y(q)
    end do
end do

do p=1,NAO
    do q=1,NAO
        H31(p,q) = conjg(H13(q,p))  ! initialize complex 0.0  H31 is the Hermitian conjuagate of H13
    end do
end do 



do p=1,NAO
    do q=1,NAO
        H22(p,q) = (0.0,0.0)  ! initialize complex 0.0  H31 is the Hermitian conjuagate of H13
    end do
end do 

do p = 1, NAO
    do q = 1, NAO
        H22(p,q) = H22(p,q) + &
            (conjg(x(p)) * conjg(x(q)) * W(p,q) * x(p) * x(q)) / 4.0 - &
            (conjg(x(p)) * conjg(y(q)) * W(p,q) * x(p) * y(q)) / 4.0 - &
            (conjg(x(q)) * conjg(y(p)) * W(p,q) * x(q) * y(p)) / 4.0 + &
            conjg(x(q)) * conjg(y(q)) * V(p,q) * x(p) * y(p) + &
            (conjg(y(p)) * conjg(y(q)) * W(p,q) * y(p) * y(q)) / 4.0
    end do
end do

do p=1,NAO
    do q=1,NAO
        H04(p,q)= (0.0,0.0)
    end do
end do


do p = 1, NAO
    do q = 1, NAO
        H04(p,q) = H04(p,q) + &
            conjg(y(p)) * conjg(y(q)) * W(p,q) * x(p) * x(q) - &
            (conjg(y(q))**2) * V(p,q) * x(p)**2
    end do
end do

do p=1,NAO
    do q=1,NAO
        H40(p,q)= conjg(H04(q,p))
    end do
end do


do p=1,NAO
    H00(p)=0.0
end do

!single sum over p
do p = 1, NAO
   H00(p) =H00(p) + (conjg(y(p))**2) * V(p,p) * (y(p)**2) + &
                     2.0 * conjg(y(p)) * h(p) * y(p)
end do
!double sum over p and q
do p = 1, NAO
    do q = 1, NAO
       H00(p) =H00(p) + conjg(x(q)) * conjg(y(q)) * V(p,q) * x(p) * y(p) + &
                          conjg(y(p)) * conjg(y(q)) * W(p,q) * y(p) * y(q)
    end do
end do

do p=1,NAO
    H02(p)=0.0
end do
 
! First sum over p
do p = 1, NAO
    H02(p) = H02(p) + 2.0 * (conjg(y(p))**2) * V(p,p) * x(p) * y(p) + &
                     2.0 * conjg(y(p)) * h(p) * x(p)
end do

! Double sum over p and q
do p = 1, NAO
    do q = 1, NAO
        H02(p) = H02(p) + conjg(x(q)) * conjg(y(q)) * V(p,q) * x(p)**2 - &
                          (conjg(y(p))**2) * V(q,p) * x(q) * y(q) + &
                          conjg(y(p)) * conjg(y(q)) * W(p,q) * x(p) * y(q) + &
                          conjg(y(p)) * conjg(y(q)) * W(q,p) * x(p) * y(q)
    end do
end do

do p=1,NAO
    H20(p) = conjg(H02(p))
end do
