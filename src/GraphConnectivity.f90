Module GraphConnMod
Use Vars, Only: MaxAt,Numat,C,RbondMax,NA

Integer(4) Marked(MaxAt)
Integer(4) iConnectivity,nmarked
Logical(1) lConn(MaxAt,MaxAt)

End module
!******************************************************
Subroutine AdjMatrix
Use Vars, Only: iAMx
Use GraphConnMod
Use Elements, Only: RA
Implicit Real(8) (A-H,O-Z)

If (rBondMax>0.d0) rbond=rBondMax
iAMx=0
nbndx=0

Do i=2,Numat
    Do j=1,i-1
        rij=Distance(i,j,Numat,C)
        If (rBondMax<0.d0) rbond=(RA(NA(i))+RA(NA(j)))*1.15d0
        If (rij<rbond) Then
            iAMx(i,j)=1
            iAMx(j,i)=1
            nbndx=nbndx+1
        Endif
    Enddo
Enddo

End
!******************************************************
Subroutine GraphConnectivity(iConnectivity0)
Use Vars, Only: iAMx,nbndx,StrDiameter
Use GraphConnMod
Use Elements, Only: RA
Implicit Real(8) (A-H,O-Z)

If (rBondMax>0.d0) rbond=rBondMax
lConn=0
iAMx=0
nbndx=0
StrDiameter=0.d0

Do i=2,Numat
    Do j=1,i-1
        rij=Distance(i,j,Numat,C)
        If (rij>StrDiameter) StrDiameter=rij
        If (rBondMax<0.d0) rbond=(RA(NA(i))+RA(NA(j)))*1.15d0
        If (rij<rbond) Then
            lConn(i,j)=.true.
            lConn(j,i)=.true.
            iAMx(i,j)=1
            iAMx(j,i)=1
            nbndx=nbndx+1
        Endif
    Enddo
Enddo

iConnectivity=0
Marked=0
nMarked=0
Do i=1,Numat
    If (Marked(i)==0) Then
        iConnectivity=iConnectivity+1
        Call MarkAtoms(i)
    Endif
    If (nMarked==Numat) Exit
Enddo

iConnectivity0=iConnectivity

End
!***************************************************************
Recursive Subroutine MarkAtoms(iHere)
Use GraphConnMod
Implicit Real(8) (A-H,O-Z)

If (Marked(iHere)>0) Return

Marked(iHere)=iConnectivity
nmarked=nmarked+1
Do iNext=1,Numat
    If (iNext==iHere) Cycle
    If (lConn(iHere,iNext).and.Marked(iNext)==0) Then
        Call MarkAtoms(iNext) 
        If (nmarked==Numat) Return
    Endif
Enddo

End