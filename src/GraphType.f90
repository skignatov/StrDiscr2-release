Subroutine SetupGraphs
Use Vars

Implicit Real(8) (A-H,O-Z)

Character(255) Str,FilGraph
Real(8) A(Numat,Numat),SignatureX(Numat),Signature(Numat),V(Numat,Numat),RA(Numat)
Integer(4) iValX(Numat),iVal(Numat),IA(numat)
Character(10) SubStr(MaxSubStr),buf10
Logical lexist

! Set up adjacency matrices (third index - Number of non-isomorphic connected graphs from McKay page http://users.cecs.anu.edu.au/~bdm/data/graphs.html)
Allocate(iAMi(Numat,Numat),iAMx(Numat,Numat))
If (Numat<10) Then
    ng=N2NG(Numat)
    Allocate(iAM(Numat,Numat,ng))
    Allocate(nEdges(ng))
    If (iGraphRecognition==0) Return
    Select Case (Numat)
    Case (1)
        ng=1
        iAMi(1,1)=0
        nEdges(1)=0
    Case(2)
        ng=1
        iAM(1:2,1,1)=(/0,1/)
        iAM(1:2,2,1)=(/1,0/)
        nEdges(1)=1
    Case(3)
        ng=2
        iAM(1:3,1,1)=(/0,1,0/)
        iAM(1:3,2,1)=(/1,0,1/)
        iAM(1:3,3,1)=(/0,1,0/)
        nEdges(1)=2
        iAM(1:3,1,2)=(/0,1,1/)
        iAM(1:3,2,2)=(/1,0,1/)
        iAM(1:3,3,2)=(/1,1,0/)
        nEdges(2)=3
    Case(4)
        ng=6
        iAM(1:4,1,1)=(/0,0,0,1/)
        iAM(1:4,2,1)=(/0,0,0,1/)
        iAM(1:4,3,1)=(/0,0,0,1/)
        iAM(1:4,4,1)=(/1,1,1,0/)
        nEdges(1)=3
        iAM(1:4,1,2)=(/0,0,1,1/)
        iAM(1:4,2,2)=(/0,0,0,1/)
        iAM(1:4,3,2)=(/1,0,0,0/)
        iAM(1:4,4,2)=(/1,1,0,0/)
        nEdges(2)=3
        iAM(1:4,1,3)=(/0,0,1,1/)
        iAM(1:4,2,3)=(/0,0,0,1/)
        iAM(1:4,3,3)=(/1,0,0,1/)
        iAM(1:4,4,3)=(/1,1,1,0/)
        nEdges(3)=4
        iAM(1:4,1,4)=(/0,0,1,1/)
        iAM(1:4,2,4)=(/0,0,1,1/)
        iAM(1:4,3,4)=(/1,1,0,0/)
        iAM(1:4,4,4)=(/1,1,0,0/)
        nEdges(4)=4
        iAM(1:4,1,5)=(/0,0,1,1/)
        iAM(1:4,2,5)=(/0,0,1,1/)
        iAM(1:4,3,5)=(/1,1,0,1/)
        iAM(1:4,4,5)=(/1,1,1,0/)
        nEdges(5)=5
        iAM(1:4,1,6)=(/0,1,1,1/)
        iAM(1:4,2,6)=(/1,0,1,1/)
        iAM(1:4,3,6)=(/1,1,0,1/)
        iAM(1:4,4,6)=(/1,1,1,0/)
        nEdges(6)=6
    Case default
    iDatFormat=2
    If (iDatFormat==0) Then
        Write(FilGraph,'(''graph'',i1,''c.txt'')')Numat
        Inquire(File=FilGraph,EXIST=lexist)
        If (.not.lexist) Then
            Write(6,'(/'' ERROR! Cannot find canonical graph patterns file '',a255)')FilGraph
            Write(6,'( ''        Provide proper file or use the keyword NOGRAPH to reject graph-type recognition.'')')
            Stop
        Endif
        Open(1,File=FilGraph)
        ig=0
        Do While (.not.EOF(1))
            Read(1,'(a255)')Str
            If (INDEX(Str,'Graph')==1) Then
                ig=ig+1
                Do i=1,Numat
                    Read(1,'(<Numat>i1)')iAMi(i,1:Numat)
                Enddo
                nEdges(ig)=Sum(iAMi)/2
                iAM(1:Numat,1:Numat,ig)=iAMi
            Endif
        Enddo
        Close(1)
    ElseIf (iDatFormat==1) Then
        Write(FilGraph,'(''graph'',i1,''.dat'')')Numat
        Inquire(File=FilGraph,EXIST=lexist)
        If (.not.lexist) Then
            Write(6,'(/'' ERROR! Cannot find canonical graph patterns file '',a255)')FilGraph
            Write(6,'( ''        Provide proper file or use the keyword NOGRAPH to reject graph-type recognition.'')')
            Stop
        Endif
        Open(1,File=FilGraph)
        ig=0
        Do While (.not.EOF(1))
            Read(1,'(a255)')Str
                    write(6,'(a255)')Str
            If (INDEX(Str,'****')>0) Then
                ig=ig+1
                Do i=1,Numat
                    Read(1,'(a255)')Str
                    write(6,'(i5,1x,a255)')ig,Str
                    If (ig==44) Then
                       x=0.d0
                    Endif
                    Call SubString(Str,MaxSubStr,nsubstr,SubStr)
                    Do j=1,nsubstr
                        buf10=AdjustL(SubStr(j))
                        buf10=AdjustR(buf10)
                        Write(6,'(a10,2i5)')buf10,i,j
                        Read(buf10,'(i10)')itmp
                        iAMi(i,i+j-1)=itmp
                        If (j>1) iAMi(i+j-1,i)=itmp
                    Enddo
                Enddo
                nEdges(ig)=Sum(iAMi)/2
                iAM(1:Numat,1:Numat,ig)=iAMi
            Endif
        Enddo
        Close(1)
        Write(6,'(/i10,'' canonical graph structures were read in from file '',a<Len_Trim(FilGraph)>)')ng,Trim(FilGraph)
    Endif
    End Select

    
Return
    
! Check for isomorphism
    n=Numat
    ii=0
    Do i=1,ng-1
        iAMi=iAM(1:n,1:n,i)
        nbndi=Sum(iAMi)/2
        Do j=1,n
            iVal(j)=Sum(iAmi(j,1:n))
        Enddo
        RA=Dble(iVal)
        Call HeapSort1(n,ra,ia)
        iVal=INT(RA)
        A=Dble(iAMi)
        Call JacobiSorted(A,n,n,Signature,V,NRrot,1)
        Do j=i+1,ng
            iAMx=iAM(1:n,1:n,j)
            nbndx=Sum(iAMi)/2
            If (nbndx/=nbndi) Cycle
            Do jj=1,n
                iValX(jj)=Sum(iAMx(jj,1:n))
            Enddo
If (Sum(iValX)==0) Then
    itmp=1
Endif
            RA=Dble(iValX)
            Call HeapSort1(n,ra,ia)
            iValX=INT(RA)
            id=Sum(IABS(iValX-iVal))
            If (id/=0) Cycle
            A=Dble(iAMx)
            Call JacobiSorted(A,n,n,SignatureX,V,NRrot,1)
            xx=Sum(DABS(Signature-SignatureX))
            If (xx<1.d-6) Then
                ii=ii+1
                Write(6,'(/'' WARNING! Canonical graphs '',i8,'' and'',i8,'' can be isomorphic! This warning appeared'',i5,'' times'')')i,j,ii
                Write(6,'('' Valences  '',i8,'' : '',<n>i15)')i,iVal
                Write(6,'('' Valences  '',i8,'' : '',<n>i15)')j,iValX
                Write(6,'('' Signature '',i8,'' : '',<n>f15.8)')i,Signature
                Write(6,'('' Signature '',i8,'' : '',<n>f15.8)')j,SignatureX
                Write(6,'('' Adjacency matrices of '',i8,'' and'',i8)')i,j
                Do k=1,n
                    Write(6,'(<n>i2,5x,<n>i2)')(iAMi(k,1:n),iAMx(k,1:n))
                Enddo
            Endif
        Enddo
    Enddo
    Write(6,*)
Endif
    
End
!*************************************************************************
Subroutine GraphString
Use Vars, Only: iAMx,g6str,Numat
Implicit Real(8) (A-H,O-Z)

Call AM2g6(Numat,iAMx,g6str)

End
!*************************************************************************
Subroutine DetermineGraphType(istr,igtype)
Use Vars

Implicit Real(8) (A-H,O-Z)

Integer(4) iA(0:Numat,0:Numat),ians(0:Numat),iAMy(Numat,Numat)
Real(8) A(Numat,Numat),SignatureX(Numat),Signature(Numat),V(Numat,Numat),RA(Numat)
Integer(4) ans(Numat)
Integer(4) iValX(Numat),iVal(Numat),IIA(numat)

igtype=-1
If (Numat>9) Return

iPrint6=0
!If (istr==141) iPrint6=1
!If (istr==144) iPrint6=1
!If (istr==145) iPrint6=1

n=Numat
m=Numat
ng=N2NG(NUmat)
iA=0
nbndx=Sum(iAMx)/2

Call AM2g6(n,iAMx,g6str)

Do i=1,n
    iValX(i)=Sum(iAMx(i,1:n))
Enddo
RA=Dble(iValX)
Call HeapSort1(n,RA,IIA)
iValX=INT(RA)

A=Dble(iAMx)
SignatureX=0.d0
Call JacobiSorted(A,n,n,SignatureX,V,NRrot,1)


xxmin=999999.d99
Do ig=1,ng
    iAMi=iAM(1:n,1:n,ig)

    ! Check for number of edges
    nbndi=Sum(iAMi)/2
    If (nbndx/=nbndi) Cycle

    ! Check for vertex valences
    Do i=1,n
        iVal(i)=Sum(iAMi(i,1:n))
    Enddo
    RA=Dble(iVal)
    Call HeapSort1(n,RA,IIA)
    iVal=INT(RA)
    idval=Sum(IABS(iVal-iValX))
    If (idval/=0) Cycle    
    
    If (iPrint6>0) Then
        Write(6,*)
        Do i=1,Numat
            Write(6,'(<n>i2,5x,<n>i2,5x,<n>i2)')iAMx(i,1:n),iAMi(i,1:n)
        Enddo
    Endif
    
    ! Check for adjacency matrix signatures
    A=Dble(iAMi)
    Signature=0.d0
    Call JacobiSorted(A,n,n,Signature,V,NRrot,1)
    xx=Sum(DABS(Signature-SignatureX))
    If (xx<xxmin) Then
        igbest=ig
        xxmin=xx
    Endif
    If (iPrint6>0) Then
        Write(6,'(/4i5,2f12.6)')istr,ig,nbndi,nbndx,xx
        Write(6,'(''Graph X:'',<n>f12.6)')SignatureX
        Write(6,'(''Graph i:'',<n>f12.6)')Signature
    Endif
   
Enddo 

If (iPrint6>0) Then
    Write(6,'(/'' Structure'',i5,''  --> graph'',i5,''  Signature discrepancy:'',f15.8)')istr,igbest,xxmin
    Do i=1,Numat
        Write(6,'(<n>i2,8x,<n>i2)')iAMx(i,1:n),iAM(i,1:n,igbest)
    Enddo
    Write(6,'(//)')
Endif

igtype=igbest


End
!*************************************************************************
Subroutine Hungarian(n,m,a0,ans0,cost)
Implicit Real(8) (A-H,O-Z)

Real(8) a0(n,m),cost
Integer(4) ans0(n)

Real(8) a(0:n,0:m),u(0:n),v(0:m),minv(0:m),delta,cur
Integer(4) p(0:m),way(0:m),ans(0:n)
Logical used(0:m)
Real(8),parameter:: XINF=999999.d99 ! large value, should be greater than any element of A0


! Improved Hungarian algorithm for solving of an Assignment Problem with scaling as O(n**3) / O(n**2 * m)  (version for Real(8) data)
! Translated to FORTRAN by S.K.Ignatov from C++ sample code of Andrey Lopatin described at https://e-maxx.ru/algo/assignment_hungary 
!   n - number of workers 
!   m - number of jobs (m>=n)
!   A0(n,m) - cost matrix
!   ans0(n) - assignments for each worker

a=0.d0
a(1:n,1:m)=a0(1:n,1:m)

u=0.d0
v=0.d0

p=0
way=0


Do i=1,n
	p(0) = i
	j0 = 0
	minv=XINF
    used=.False.
	Do 
		used(j0) = .true.
		i0 = p(j0)
        delta = XINF
		Do j=1,m
			If (.not.used(j)) Then
				cur = a(i0,j)-u(i0)-v(j)
				If (cur < minv(j)) Then
					minv(j) = cur
                    way(j) = j0
                Endif
				If (minv(j) < delta) Then
					delta = minv(j)
                    j1 = j
                Endif
            EndIf
        Enddo
		Do j=0,m
			If (used(j)) Then
				u(p(j)) = u(p(j)) + delta
                v(j) = v(j) - delta
			Else
				minv(j) = minv(j) - delta
            Endif
        Enddo
		j0 = j1
        If (.not.(p(j0)/=0)) Exit
    Enddo
	Do 
		j1 = way(j0)
		p(j0) = p(j1)
		j0 = j1
        If (.not.(j0>0)) Exit
    Enddo
Enddo
    
Do j=1,m
	ans(p(j)) = j
Enddo

Cost=-v(0)
ans0(1:n)=ans(1:n)

End
!*************************************************
Subroutine HungarianInt(n,m,a,ans,cost)
Implicit Real(8) (A-H,O-Z)

Integer u(0:n), v(0:m), p(0:m), way(0:m)
Integer(4) minv(0:m)
Logical used(0:m)
Integer(4) delta,cur
Integer(4) a(0:n,0:m)
Integer(4) ans(0:m),cost

INF=100
u=0
v=0
p=0
way=0
ans=0

Do i=1,n
	p(0) = i
	j0 = 0
	minv=INF
    used=.False.
	Do 
		used(j0) = .true.
		i0 = p(j0)
        delta = INF     !,  j1;
		Do j=1,m
			If (.not.used(j)) Then
				cur = a(i0,j)-u(i0)-v(j)
				If (cur < minv(j)) Then
					minv(j) = cur
                    way(j) = j0
                Endif
				If (minv(j) < delta) Then
					delta = minv(j)
                    j1 = j
                Endif
            EndIf
        Enddo
		Do j=0,m
			If (used(j)) Then
				u(p(j)) = u(p(j)) + delta
                v(j) = v(j) - delta
			Else
				minv(j) = minv(j) - delta
            Endif
        Enddo
		j0 = j1
        If (.not.(p(j0)/=0)) Exit
    Enddo
!	} while (p[j0] != 0);
	Do 
		j1 = way(j0)
		p(j0) = p(j1)
		j0 = j1
!	} while (j0);
        If (.not.(j0>0)) Exit
    Enddo
Enddo
    
Do j=1,m
	ans(p(j)) = j
Enddo

Cost=-v(0)

    End
!************************************************************************
Subroutine AM2g6(n,iAM,g6str)
Implicit Real(8) (A-H,O-Z)

Integer(4) iAM(n,n)
Character(1000) g6str
Integer(4) iX(6,1000)
Character(6) buf

k=0
kg=1
Do j=2,n
    Do i=1,j-1
        k=k+1
        If (k>6) Then
            k=1
            kg=kg+1
        Endif
        iX(k,kg)=iAM(i,j)
    Enddo
Enddo
If (k<6) iX(k+1:6,kg)=0

g6str(1:1)=CHAR(n+63)
Do i=1,kg
    Write(buf,'(6i1)')iX(1:6,i)
    Read(buf,'(b6)')itmp
    itmp=itmp+63
    g6str(i+1:i+1)=CHAR(itmp)
Enddo

End
!************************************************************************
Subroutine g62am(g6str,n,iAM)
Implicit Real(8) (A-H,O-Z)

Integer(4) iAM(n,n)
Character(1000) g6str
Integer(4) iX(1000)
Character buf*6,ch*1

!n=ICHAR(g6str(1:1))-63
ll=Len_Trim(g6str)
k=0
Do i=2,ll
    ch=g6str(i:i)
    ich=ICHAR(ch)-63
    Write(buf,'(b6)')ich
    Do j=1,6
        ch=buf(j:j)
        Read(ch,'(i1)')itmp
        k=k+1
        iX(k)=itmp
    Enddo
Enddo

iAM=0
k=0
Do j=2,n
    Do i=1,j-1
        k=k+1
        iAM(i,j)=iX(k)
        iAM(j,i)=iX(k)
    Enddo
Enddo

End