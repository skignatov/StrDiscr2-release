Subroutine SymGroup(kRotTyp,iSym,iSymGrp,NRotSym,SymGrp)
USE Vars, Only: N=>Numat,NUC=>NA,CC0=>C,AMass0=>AMass,SymEps,Angs

Implicit Real(8) (A-H,O-Z)

Real(8) C0(3,N),C(3,N),AMass(N)
Integer(4) NA(N)

Real(8) CM(3),PMOI(3),PAxes(3,3),T(3,3)
Integer(4) LowSymEl(7)
Character(4) SymGrp

! Subroutine SymGroup determines the symmetry group of molecule
! defined by their coordinates C0, and masses AMass


SymGrp='    '
iSymGrp=0
Do i=1,N
	Do k=1,3
		C0(k,i)=CC0(k,i)*Angs
	Enddo
	NA(i)=NUC(i)
	AMass(i)=AMass0(i)
Enddo
!IF (SymEps/=0.005d0) Write(6,'(//'' WARNING! Non-standard value of SymEps is used ''/ &
!								 '' to recognize the molecular symmetry:'',f20.10)')SymEps

!
! Reduce the molecule to principal axes
!

!CALL Inertia(N,C0,AMass,CM,PMOI,PAxes,kRotTyp)
Call InertiaNew(2,-1,1,N,C0,AMass,CM,PMOI,PAxes,kRotTyp)
!
! Determine the Symmetry operations
!
1 Continue
isym=kRotTyp*10**7
If (kRotTyp==5) Then		! Case of Linear molecule
	Do k=1,3
		Call SymOp1(k,T)
		Call Equiv1(N,NA,C0,T,lse)
		LowSymEl(k)=lse
	Enddo
	Call PermuteAxes(N,C0,LowSymEl)
	Do k=1,7
		Call SymOp1(k,T)
		Call Equiv1(N,NA,C0,T,lse)
		LowSymEl(k)=lse
	Enddo

ElseIf (kRotTyp==4) Then	! Case of non-linear Low symmetry 
	! Check for rotation symmetry
	Do k=1,3
		Call SymOp1(k,T)
		Call Equiv1(N,NA,C0,T,lse)
		LowSymEl(k)=lse
	Enddo
	Call PermuteAxes(N,C0,LowSymEl)
	! Check for all the low symmetry elements
	Do k=1,7
		Call SymOp1(k,T)
		Call Equiv1(N,NA,C0,T,lse)
		LowSymEl(k)=lse
!		Write(6,*)k,lse
	Enddo

ElseIf (kRotTyp==2.or.kRotTyp==3) Then		! Medium symmetry
	LowSymEl=0
	! Permute axes for oblate top
	If (kRotTyp==3) Then
		LowSymEl(1)=1
		Call PermuteAxes(N,C0,LowSymEl)
		LowSymEl(1)=0
	Endif
	! Check for Cn
	Do i=2,8
		Call SymOp2(1,i,T)
		Call Equiv1(N,NA,C0,T,lse)
		If (lse==1) Then
			LowSymEl(1)=i
			LowSymEl(2)=1
		Endif
	Enddo
	If (LowSymEl(1)==0) Then				! Occasional degeneracy
		!Write(6,'('' Occasional degeneracy'')')
		kRotTyp=4
		goto 1
	Endif
	! Check for S2n
	Call SymOp2(2,2*LowSymEl(1),T)
	Call Equiv1(N,NA,C0,T,LowSymEl(3))
	! Check for SigmaH
	Call SymOp2(3,0,T)
	Call Equiv1(N,NA,C0,T,LowSymEl(4))
	! Check for SigmaV
	Do i=0,N
		If (i==0) Then
			C=C0
			Goto 20
		Endif
		RR=DSQRT(C0(1,i)**2+C0(2,i)**2)
		If (RR<0.01d0) Cycle
		ct=C0(1,i)/RR
		st=C0(2,i)/RR
		Call RotXYZ(3,-1,-1,ct,st,N,C0,C)
	20  Continue
		Call SymOp2(4,0,T)
		Call Equiv1(N,NA,C,T,lse)
		If (lse==1) Then
			LowSymEl(5)=1
			C0=C
			Exit
		Endif
	Enddo
ElseIf (kRotTyp==1) Then
	
	Call HighSymmetry(kRotTyp,N,NA,C0,LowSymEl,iDegen)
	If (iDegen>0) Then
		Write(6,'('' Occasional degeneracy'')')
		kRotTyp=2
		Goto 1
	Endif

Endif

! Recognize the symmetry group
!Write(6,'('' PMOI: '',3f10.6)')(PMOI(k),k=1,3)
Do k=1,7
	iSym=isym+LowSymEl(k)*10**(7-k)
Enddo
If (kRotTyp==4.or.kRotTyp==5) Then
	If     (iSym==40000000) Then ! C1
		iSymGrp=1
		SymGrp='C1  '
		NRotSym=1
	ElseIf (iSym==40000001) Then ! Ci
		iSymGrp=2
		SymGrp='Ci  '
		NRotSym=1
	ElseIf (iSym==40000010.or.iSym==40000100.or.iSym==40001000) Then ! Cs
		iSymGrp=3
		SymGrp='Cs  '
		NRotSym=1
	ElseIf (iSym==40010000) Then ! C2
		iSymGrp=4
		SymGrp='C2  '
		NRotSym=2
	ElseIf (iSym==40011100) Then	! C2v
		iSymGrp=5
		SymGrp='C2v '
		NRotSym=2
	ElseIf (iSym==40010010.or.iSym==40010011) Then ! C2h
		iSymGrp=6
		SymGrp='C2h '
		NRotSym=2
	ElseIf (iSym==41111111) Then ! D2h
		iSymGrp=7
		SymGrp='D2h '
		NRotSym=4
	ElseIf (iSym==50011100) Then ! C*v
		iSymGrp=8
		SymGrp='C*v '
		NRotSym=1
	ElseIf (iSym==51111111) Then ! D*h
		iSymGrp=9
		SymGrp='D*h '
		NRotSym=2
	Endif
ElseIf (kRotTyp==2.or.kRotTyp==3) Then
	iSym2=0
	Do k=2,5
		iSym2=iSym2+LowSymEl(k)*10**(5-k)
	Enddo
	If     (iSym2==1000) Then ! Cn
		iSymGrp=10
		Write(SymGrp,'(''C'',i1,''  '')')LowSymEl(1)
		NRotSym=LowSymEl(1)
	ElseIf (iSym2==1001) Then ! Cnv
		iSymGrp=11
		Write(SymGrp,'(''C'',i1,''v '')')LowSymEl(1)
		NRotSym=LowSymEl(1)
	ElseIf (iSym2==1010) Then ! Cnh=Dn
		iSymGrp=12
		Write(SymGrp,'(''C'',i1,''h '')')LowSymEl(1)
		NRotSym=LowSymEl(1)
	ElseIf (iSym2==1011) Then ! Dnh
		iSymGrp=13
		Write(SymGrp,'(''D'',i1,''h '')')LowSymEl(1)
		NRotSym=LowSymEl(1)*2
	ElseIf (iSym2==1101) Then ! Dnd
		iSymGrp=14
		Write(SymGrp,'(''D'',i1,''d '')')LowSymEl(1)
		NRotSym=LowSymEl(1)*2
	Endif
	If (iSymGrp==0) Then
		kRotTyp=4
		Goto 1
	Endif
ElseIf (kRotTyp==1) Then
	iSym3=0
	Do k=1,5
		iSym3=iSym3+LowSymEl(k)*10**(5-k)
	Enddo
	If     (iSym3==11000) Then ! T
		iSymGrp=15
		SymGrp='T   '
		NRotSym=12	
	ElseIf (iSym3==11001) Then ! Td
		iSymGrp=16
		SymGrp='Td  '
		NRotSym=12
	ElseIf (iSym3==11012) Then ! Th
		iSymGrp=17
		SymGrp='Th  '
		NRotSym=12	
	ElseIf (iSym3==11102) Then ! O
		iSymGrp=18
		SymGrp='O   '
		NRotSym=24
	ElseIf (iSym3==11112) Then ! Oh
		iSymGrp=19
		SymGrp='Oh  '
		NRotSym=24	
		Write(6,'(''Oh'')')
	ElseIf (iSym3==10222) Then ! I
		iSymGrp=20
		SymGrp='Ih  '	! Really, this is group I
		NRotSym=60	
	Endif
	If (iSymGrp==0) Then
		kRotTyp=2
		Goto 1
	Endif
ElseIf (kRotTyp==6) Then
	iSym=60000000
	iSymGrp=21
	SymGrp='O(3)'
	NRotSym=1
Endif
		

End
!**************************************************************
Subroutine HighSymmetry(kRotTyp,N,NA,C0,LowSymEl,iDegen)
Implicit Real(8) (A-H,O-Z)

Real(8) C0(3,N),C(3,N),C1(3,N),T(3,3)
Integer(4) NA(N),LowSymEl(7)
Integer(4) IEqAt(N,0:N),iAtC3(N)

! Subroutine HighSymmetry recognizes the high symmetry groups: 
! T, Td, Th, O, Oh, I
! If iDegen>0, the occasional degeneracy is occured.

Eps=0.005d0
Pi=3.141592653589793d0		
ta=109.4712d0*Pi/180.d0
iDegen=0

! First, equivalent atoms are sought
IEqAt=0
mink=0
Do i=1,N
	nai=NA(i)
	Ri=C0(1,i)**2+C0(2,i)**2+C0(3,i)**2
	k=0
	Do j=1,N
		naj=NA(j)
		If (naj/=nai) Cycle
		Rj=C0(1,j)**2+C0(2,j)**2+C0(3,j)**2
		If (DABS(Ri-Rj)>Eps) Cycle
		k=k+1
		IEqAt(i,k)=j
	Enddo
	iEqAt(i,0)=k
	If (mink==0.and.k>=4) mink=i 
	If (k>=4.and.k<mink) mink=i
Enddo
!Do i=1,N
!	Write(6,'(10i3)')(iEqAt(i,j),j=0,N)
!Enddo
If (mink==0) Then
	iDegen=1
	Return
Endif

! Search for C3 elements
km=IEqAt(mink,0)
lC3=0
Do i=1,km
	ia=IEqAt(mink,i)
	Call RotateA2Z(ia,N,C0,C)
	Call SymOp2(1,3,T)
	Call Equiv1(N,NA,C,T,lse)
	If (lse==1) Then
		lC3=lC3+1
!		LowSymEl(1)=lC3
		iAtC3(lC3)=ia
	Endif
Enddo
If (lC3<4) Then
	iDegen=2
	Return
EndIf
LowSymEl(1)=1

! Check for tetrahedral angle between C3 axes
itd=0
ita=0
itb=0
Do i=1,lC3
	ia=iAtC3(i)
	a2=C0(1,ia)**2+C0(2,ia)**2+C0(3,ia)**2
	Do j=i+1,lC3
		ja=iAtC3(j)
		b2=C0(1,ja)**2+C0(2,ja)**2+C0(3,ja)**2
		c2=(C0(1,ia)-C0(1,ja))**2+(C0(2,ia)-C0(2,ja))**2+(C0(3,ia)-C0(3,ja))**2
		ca=0.5d0*(a2+b2-c2)/DSQRT(a2*b2)
		Ang=DACOS(ca)
		If (DABS(Ang-ta)<0.01d0) Then
			itd=1
			If (ita==0) Then
				ita=ia
				itb=ja
			Endif
		Endif
	Enddo
Enddo
If (itd==0) Then		! Icosahedral symmetry
	LowSymEl(2)=0
	LowSymEl(3:5)=2
	Return
Endif
LowSymEl(2)=1

! Orient C3 axes for cubic groups
ia=ita
ib=itb
CALL RotateA2Z(ia,N,C0,C)	
Rxy=DSQRT(C(1,ib)**2+C(2,ib)**2)
ct=C(1,ib)/Rxy
st=C(2,ib)/Rxy
Call RotXYZ(3,-1,-1,ct,st,N,C,C1)
Ang=Pi/2.d0-0.5d0*ta
Call RotXYZ(2,1,1,Ang,Ang,N,C1,C)
Call RotXYZ(1,1,0,45.d0,0.d0,N,C,C1)
!Call PrintC(N,NA,C1)

! Check for C2d axes
Call RotXYZ(3,1,0,45.d0,0.d0,N,C1,C)
Call SymOp1(1,T)
Call Equiv1(N,NA,C,T,LowSymEl(3))
! Check for SigmaD
Call SymOp1(4,T)
Call Equiv1(N,NA,C,T,LowSymEl(5))
! Check for i center
Call SymOp1(7,T)
Call Equiv1(N,NA,C,T,LowSymEl(4))
Call RotXYZ(3,-1,0,45.d0,0.d0,N,C,C0)

If (LowSymEl(3)==1.or.LowSymEl(4)==1) LowSymEl(5)=2


End
!**************************************************************
Subroutine RotateA2Z(iAtom,N,C0,C)
Implicit Real(8) (A-H,O-Z)

Real(8) C0(3,N),C(3,N),T(3,3),C1(3,N)

x=C0(1,iAtom)
y=C0(2,iAtom)
z=C0(3,iAtom)
Rxy=DSQRT(x*x+y*y)
R=DSQRT(Rxy*Rxy+z*z)
If (Rxy>1.d-8) Then
	cv=x/Rxy
	sv=y/Rxy
Else
	cv=1.d0
	sv=0.d0
Endif
If (R>1.d-8) Then
	ct=z/R
	st=DSQRT(1.d0-ct*ct)
Else
	ct=1.d0
	st=0.d0
Endif

Call RotXYZ(3,-1,-1,cv,sv,N,C0,C1)
Call RotXYZ(2,-1,-1,ct,st,N,C1,C)

End
!**************************************************************
Subroutine RotXYZ(iAxis,iDir,iRad,A,B,N,C0,C)
Implicit Real(8) (A-H,O-Z)

Real(8) C0(3,N),C(3,N),T(3,3)

! Subroutine RotXYZ rotates coordinates C0 around axis iAxis (OX, OY, or OZ)
! by the defined angle A (Degrees if iRad==0, radians if iRad/=0)
! or by cosinus and sinus A and B of the angle if iRad<0
! The rotation is a subject of right-hand rule if iDir>=0, and back otherwise.

Pi=3.141592653589793d0		

If (iRad>=0) Then
	X=A
	If (iRad==0) X=A*Pi/180.d0
	ct=DCOS(X)
	st=DSIN(X)
Else
	ct=A
	st=B
Endif

If (iDir<0) st=-st

T=0.d0
If (iAxis==3) Then
	T(1,1)=ct
	T(2,2)=ct
	T(3,3)=1.d0
	T(1,2)=-st
	T(2,1)=-T(1,2)
ElseIf (iAxis==1) Then
	T(1,1)=1.d0
	T(2,2)=ct
	T(3,3)=ct
	T(2,3)=-st
	T(3,2)=-T(2,3)
	T(1,1)=1.d0
ElseIf (iAxis==2) Then
	T(1,1)=ct
	T(2,2)=1.d0
	T(3,3)=ct
	T(1,3)=+st
	T(3,1)=-T(1,3)
Endif

Do i=1,N
	Do k=1,3
		C(k,i)=T(k,1)*C0(1,i)+T(k,2)*C0(2,i)+T(k,3)*C0(3,i)
	Enddo
Enddo

End
!**************************************************************
Subroutine SymOp1(k,T)
Implicit Real(8) (A-H,O-Z)

Real(8) T(3,3) 

T=0.d0
ForAll(i=1:3) T(i,i)=-1.d0
If (k==7) Return
If (k<4) Then
	T(k,k)=1.d0
	Return
Endif
T(k-3,k-3)=1.d0
ForAll(i=1:3) T(i,i)=-T(i,i)

End
!**************************************************************
Subroutine SymOp2(k,n,T)
Implicit Real(8) (A-H,O-Z)

Real(8) T(3,3) 

T=0.d0

If (k<=2) Then				! Cn and Sn
	Pi=3.141592653589793d0		
	X=2.d0*Pi/Dble(n)
	cx=DCOS(X)
	sx=DSIN(X)
	T(1,1)=cx
	T(1,2)=-sx	
	T(2,1)=sx
	T(2,2)=cx
	If (k==1) T(3,3)=1.d0
	If (k==2) T(3,3)=-1.d0
	Return
Endif

T(1,1)=1.d0
If (k==3) Then				! SigmaH
	T(2,2)=1.d0
	T(3,3)=-1.d0
ElseIf (k==4) Then				! SigmaVx
	T(2,2)=-1.d0
	T(3,3)=1.d0
Endif

End
!**************************************************************
Subroutine Equiv1(N,NA,C,T,IsEquiv)
USE Vars, Only: SymEps
Implicit Real(8) (A-H,O-Z)


Real(8) C(3,N),T(3,3),C1(3)
Integer(4) NA(N),iPerm(N)

Eps=SymEps
IsEquiv=0
iPerm=0

IPrint=0

If (IPrint>0) Then
	Write(6,'(3f10.4)')T
	Write(6,*)
	Do i=1,N
		Write(6,'(i5,3f10.4)')NA(i),(C(j,i),j=1,3)
	Enddo
	Write(6,*)
EndIf

iEquiv=0
Do i=1,N
	Do j=1,3
		C1(j)=0.d0
		Do k=1,3
			C1(j)=C1(j)+T(j,k)*C(k,i)
		Enddo
	Enddo
	If (IPrint>0) Write(6,'(i5,3f10.4)')NA(i),(C1(j),j=1,3)
	jEquiv=0
	Do j=1,N
		If (NA(j)/=NA(i)) Cycle
		dx=DABS(C(1,j)-C1(1))
		dy=DABS(C(2,j)-C1(2))
		dz=DABS(C(3,j)-C1(3))
		If (dx>Eps.or.dy>Eps.or.dz>Eps) Cycle
		jEquiv=jEquiv+1
		iPerm(i)=j
	Enddo
	If (jEquiv>1) Then
		Write(*,'('' ERROR! Some atoms coincide after symmetry transformation!'',2i4)')i,j
		Write(*,901)
		Pause
		Return
	Endif
	iEquiv=iEquiv+jEquiv
Enddo

If (iEquiv<0.or.iEquiv>N) Then
	Write(*,'('' ERROR in Equiv1!'',i5)')iEquiv
	Write(*,901)
	Pause
	IsEquiv=0
	Return
Endif
If (iEquiv==N) IsEquiv=1
If (IPrint>0) Write(6,*)IsEquiv

901 Format(' This is probably due to the incorrect molecular geometry.'/' Attention! -- Symmetry can be incorrectly identified!.')

End
!************************************************************
Subroutine PermuteAxes(N,C0,LowSymEl)
Implicit Real(8) (A-H,O-Z)

Real(8) C0(3,N)
Integer(4) LowSymEl(7)

! Permute the axes to provide that OZ is a principal axis

	If (LowSymEl(1)==1.and.LowSymEl(2)==0.and.LowSymEl(3)==0) Then
		Do i=1,N
			tmp=C0(1,i)
			C0(1,i)=C0(2,i)
			C0(2,i)=C0(3,i)
			C0(3,i)=tmp
		Enddo
	ElseIf (LowSymEl(1)==0.and.LowSymEl(2)==1.and.LowSymEl(3)==0) Then
		Do i=1,N
			tmp=C0(3,i)
			C0(3,i)=C0(2,i)
			C0(2,i)=C0(1,i)
			C0(1,i)=tmp
		Enddo
	Endif

End