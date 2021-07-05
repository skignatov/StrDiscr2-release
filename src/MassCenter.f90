Subroutine MassCenter(mode,N,C,AM,CM,TotMass)
Implicit Real(8) (A-H,O-Z)

Real(8) C(3,N),AM(N),CM(3)

! Subroutine MassCenter calulates center of masses CM(3) and total mass of N-atom molecule
! with the atom coordinates C(3,N) and atomic masses AM(N)
!
! MODE==0 - the coordinates C(3,N) are not changed
! MODE/=0 - the coordinates C(3,N) are reduced to center of mass

TotMass=0.d0
CM=0.d0
Do i=1,N
	TotMass=TotMass+AM(i)
	Do k=1,3
		CM(k)=CM(k)+AM(i)*C(k,i)
	Enddo
Enddo
CM=CM/TotMass

If (mode==0) Return

Do i=1,N
	Do k=1,3
		C(k,i)=C(k,i)-CM(k)
	Enddo
Enddo
	
End
