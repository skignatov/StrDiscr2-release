Subroutine ReadFreqs(iu,istr,iopt,ierr,ierr2,ngeofound,netotfound,nfreq,nneg,Etot,C0,Freq,Vibr1)
Use Vars, FreqStr=>Freq, Vibr1str=>Vibr1
!Use StringMod
Implicit Real(8) (A-H,O-Z)

Character(255) Str,SubStr(MaxSubStr),Str1 !,FilLogX
Character(20) buf
Integer(4) NA0(MaxAt)
Real(8) C0(3,MaxAt),V(3*MaxAt,3*MaxAt),Freq(3*MaxAt),Vibr1(3,MaxAt),FrcConst(3*MaxAt),RedMass(3*MaxAt),XIRint(3*MaxAt)
Real(8) afreq(3),arm(3),afc(3),axir(3)

!iu=3
!Open(iu,File=FilLogX)
ierr=-1
ierr2=0
V=0.d0
iopt=0
nfreq=0
nneg=0
netotfound=0
ngeofound=0
nrm=0
nfc=0
nir=0
Do While(.not.EOF(iu))
    Read(iu,'(a255)')Str
    If (INDEX(Str,'SCF Done:')>0) Then
        Call SubString(Str,MaxSubStr,nSubStr,SubStr)
        buf=AdjustL(SubStr(5))
        Read(buf,'(f20.8)')etot
        netotfound=netotfound+1
        Cycle
    Endif
    If (INDEX(Str,'-- Stationary point found.')>0) Then
        iopt=2
        Cycle
    Endif
    If (INDEX(Str,'-- Number of steps exceeded,')>0) Then
        ierr=2
        Cycle
    Endif
    If (INDEX(Str,'Error termination via Lnk1e in')>0) Then
        If (ierr/=2) ierr=3
        Cycle
    Endif
    If (INDEX(Str,'Standard orientation:')>0) Then
!        Call ReadNXYZ(-iu,0,4,10,2,MaxAt,NA0,C0,NumAt0,AName,ierr1,Str)
        ierr1=0
        Call ReadNXYZ(-iu,0,4,10,2,MaxAt,NA0,C0,NumAt0,AName,Str)
        If (ierr1>0) Then
            ierr2=1
            Cycle
        Endif
        Numat=Numat0
        NA(1:Numat)=NA0(1:Numat)
        C(1:3,1:NUmat)=C0(1:3,1:Numat)
        ngeofound=ngeofound+1
        Cycle
    Endif
    If (Index(Str,'Frequencies --')>0) Then
        Call SubString(Str,MaxSubStr,nSubStr,SubStr)
        nfreq1=nfreq+1
        nf=nSubStr-2
        Do ii=3,nSubStr
            buf=AdjustL(SubStr(ii))
            Read(buf,'(f20.6)')ff
            nfreq=nfreq+1
            Freq(nfreq)=ff
            If (ff<0.d0) nneg=nneg+1
        Enddo
        ! Read reduced masses
        Read(iu,'(a255)')Str1
        Call SubString(Str1,MaxSubStr,nSubStr,SubStr)
        Do ii=4,nSubStr
            buf=AdjustL(SubStr(ii))
            Read(buf,'(f20.6)')rm
            nrm=nrm+1
            RedMass(nrm)=rm
        Enddo
        ! Read force constants
        Read(iu,'(a255)')Str1
        Call SubString(Str1,MaxSubStr,nSubStr,SubStr)
        Do ii=4,nSubStr
            buf=AdjustL(SubStr(ii))
            Read(buf,'(f20.6)')fc
            nfc=nfc+1
            FrcConst(nfc)=fc
        Enddo
        ! Read IR intensities
        Read(iu,'(a255)')Str1
        Call SubString(Str1,MaxSubStr,nSubStr,SubStr)
        Do ii=4,nSubStr
            buf=AdjustL(SubStr(ii))
            Read(buf,'(f20.6)')xir
            nir=nir+1
            XIRint(nir)=xir
        Enddo
        Read(iu,*)
        Do ia=1,Numat
            Read(iu,'(a255)')Str
            Call SubString(Str,MaxSubStr,nSubStr,SubStr)
            ij=2
            Do i=nfreq1,nfreq
                Do j=1,3
                    ij=ij+1
                    kk=(ia-1)*3+j
                    buf=AdjustL(SubStr(ij))
                    Read(buf,'(f20.6)')V(kk,i)
                Enddo
            Enddo
        Enddo
        ierr=0
    Endif
Enddo

kk=0
Do i=1,Numat
    Do k=1,3
        kk=kk+1
        Vibr1(k,i)=V(kk,1)
    Enddo
Enddo
            
! Printing RedMass,FrcConst,XIRint
afreq(1)=Sum(Freq(1:3))/3.d0
afreq(2)=Sum(Freq(1:6))/6.d0
afreq(3)=Sum(Freq(1:9))/9.d0

arm(1)=Sum(RedMass(1:3))/3.d0
arm(2)=Sum(RedMass(1:6))/6.d0
arm(3)=Sum(RedMass(1:9))/9.d0

afc(1)=Sum(FrcConst(1:3))/3.d0*1.d3   ! mDyn/A -> N/cm  !*143.836d0    ! Convert mDyn/A -> kcal/mol/A**2
afc(2)=Sum(FrcConst(1:6))/6.d0*1.d3
afc(3)=Sum(FrcConst(1:9))/9.d0*1.d3

Write(11,'(''!Averaged (3/6/9)values:     Freqs                           Red.Mass                        Force constants'')')
Write(11,'(''!geo '',i5,''  Numat:'',i4,2x,3(3f10.4,2x))')istr+1,Numat,afreq(1:3),arm(1:3),afc(1:3)
!Do i=1,nfreq
!    Write(11,'(i5,4f10.2)')i,Freq(i),RedMass(i),FrcConst(i),XIRint(i)
!Enddo

End
