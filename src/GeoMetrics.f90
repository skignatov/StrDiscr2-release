!************************************************************************************
Subroutine GeoMetrics(mode,Numat,NA,C,Score,D0)

Use Vars, Only: MaxAt   !,Amass
Use Elements, Only: RA

Implicit Real(8) (A-H,O-Z)

Integer(4) NA(Numat)
Real(8) C(3,Numat),S(4*Numat,4*Numat),D(4*Numat),V(4*Numat,4*Numat),CM(3)
Real(8) Alp(Numat),D0(4*Numat)
!Static Alp,D0

n4=Numat*4

! Set Gaussian exponents
!    xn=Dble(Numat)
!    bs=((Product(BoxSize))**(1.d0/3.d0))    ! Geometric mean of a,b,c
!    !Alp(1:120)=1.d0/RA(1:120)**2                       ! Set up as in Goedecker2013
!    !Alp=(Product(BoxSize)/Dble(Numat))**(-1.d0/3.d0)   !good choice
!    !Alp=2.d0/(bs*bs)
!    rb=0.5d0*bs/xn**(1.d0/3.d0)      ! Average bond radius of N atoms in Vol (average closest distance ra=(Vol/n)^(-1/3), rb=0.5*ra)
!    Alp=1/rb**2

If (mode==0) Then
    !Alp(1:Numat)=1.d0/RA(NA(1:Numat))**2
    Do i=1,Numat
        Alp(i)=1.d0/RA(NA(i))**2
    Enddo
Endif

!Call MassCenter(1,Numat,C,AMass,CM,TotMass)
Call Sint(Numat,Alp,C,S)
Call JacobiSorted(S,n4,n4,D,V,nrot,-1)

If (mode==0) Then
    D0(1:n4)=D(1:n4)                ! Remember metrics of initial structure
    Score=0.d0
    Return
Else
    D(1:n4)=D(1:n4)-D0(1:n4)
    Score=dot_product(D,D)/Dble(n4)
!    Score=DSQRT(Score)
!    Score=Score*bs
    Score=Score*1.d4
Endif

End
!************************************************************************************
Subroutine Sint(Numat,Alp,C,S)
!Use Elements, Only: RA

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,Numat),S(4*Numat,4*Numat),Alp(Numat)
Real(8), parameter:: A2au=0.529177d0

n4=Numat*4
S=0.d0

Do i=2,Numat
!    nai=NA(i)
!    ai=Alp(nai)
    ai=Alp(i)
    i0=(i-1)*4
    Do j=1,i-1
!        naj=NA(j)
!        aj=Alp(naj)
        aj=Alp(j)
        a1=2.d0*DSQRT(ai*aj)/(ai+aj)
        a2=(ai*aj)/(ai+aj)
        rij=Distance(i,j,Numat,C)
        ss=DSQRT(a1*a1*a1)*DEXP(-a2*rij*rij)
        a1ss=a1*ss
        spx=-a1ss*DSQRT(aj)*(C(1,i)-C(1,j))
        spy=-a1ss*DSQRT(aj)*(C(2,i)-C(2,j))
        spz=-a1ss*DSQRT(aj)*(C(3,i)-C(3,j))
        pxpx=a1ss*(1.d0-2.d0*a2*(C(1,i)-C(1,j))**2)
        pypy=a1ss*(1.d0-2.d0*a2*(C(2,i)-C(2,j))**2)
        pzpz=a1ss*(1.d0-2.d0*a2*(C(3,i)-C(3,j))**2)
        pxpy=-a1ss*2.d0*a2*(C(1,i)-C(1,j))*(C(2,i)-C(2,j))
        pxpz=-a1ss*2.d0*a2*(C(1,i)-C(1,j))*(C(3,i)-C(3,j))
        pypz=-a1ss*2.d0*a2*(C(2,i)-C(2,j))*(C(3,i)-C(3,j))
        j0=(j-1)*4
!
        S(i0+1,j0+1)=ss
        S(i0+1,j0+2)=spx
        S(i0+1,j0+3)=spy
        S(i0+1,j0+4)=spz
! - minus?
        S(i0+2,j0+1)=-spx
        S(i0+3,j0+1)=-spy
        S(i0+4,j0+1)=-spz
!
        S(i0+2,j0+2)=pxpx
        S(i0+3,j0+3)=pypy
        S(i0+4,j0+4)=pzpz
!
        S(i0+2,j0+3)=pxpy
        S(i0+3,j0+2)=pxpy
!
        S(i0+2,j0+4)=pxpz
        S(i0+4,j0+2)=pxpz
!
        S(i0+3,j0+4)=pypz
        S(i0+4,j0+3)=pypz
    Enddo
Enddo

Do i=1,n4
    S(i,i)=1.d0
    Do j=i,n4
        S(i,j)=S(j,i)
    Enddo
Enddo
        
End
