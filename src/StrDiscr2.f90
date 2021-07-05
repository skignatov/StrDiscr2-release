Module Vars

    Implicit Real(8) (A-H,O-Z)
    
    Integer(4), parameter::MaxAt=100,MaxStr=10000,MaxPot=1,MaxSubStr=20
    
    Real(8) C(3,MaxAt),Cstr(3,MaxAt,MaxStr),C1(3,MaxAt),C2(3,MaxAt),DM(MaxAt,MaxAt),Amass(MaxAt),PSE(MaxPot,MaxStr),Dlt(MaxPot,MaxStr)
    Integer(4) NA(MaxAt),NA1(MaxAt),Numat,NAstr(MaxAt,MaxStr),NumatStr(MaxStr),IAstr(MaxStr),Itype(0:MaxStr,MaxStr,MaxPot),NTl(MaxPot)
    Character(10) Aname(MaxAt),Aname10
    
    Real(8) epsmash,RbondMax/-1.d0/,SymEps/0.1d0/,angs/1.d0/

    Integer(1),allocatable::iAM(:,:,:)
    Integer(4),allocatable::iAMi(:,:),nEdges(:),iAMx(:,:)
    Integer(4) N2NG(11)/1,1,2,6,21,112,853,11117,261080,11716571,1006700565/    ! Number of non-isomorphic connected graphs from McKay page http://users.cecs.anu.edu.au/~bdm/data/graphs.html)
    Integer(4) nbndx
    Integer(4) IsoGraph(MaxStr),iGraphRecognition/1/
!    Character(10) Gname(:)
    
    Real(8) Freq(3,MaxAt,MaxStr),Vibr1(3,MaxAt,MaxStr),StrDiameter
    Integer(4) NFreqStr(MaxStr),NVibrStr(MaxStr)
    
    Character(1000) g6str
    
End module
!*******************************************************************
Program StrDiscr2

Use Vars
Use Elements, Only: AMS

Implicit Real(8) (A-H,O-Z)

Character(255) FilOut,Str,SearchStr,SubStr(50),FilDat,FilCan
Real(8) E0(MaxPot),E(MaxPot),dF(MaxPot),Gdx(MaxPot),X(MaxStr),CM(3),Etot(3,MaxStr),PMOI(3),PAxes(3,3),PM(3,MaxStr),dreor(0:7)
Real(8) Eaver(3,MaxStr),Estdv(3,MaxStr)

Character(255) StrLabelX(MaxStr)
Character(4) SymGrpStr(MaxStr),SymGrp

Real(8),allocatable::E1(:),E2(:)
Integer(4),allocatable::IE1(:),ITT1(:)

Real(8),allocatable::DC(:,:),XDC(:)
Real(8),allocatable::D1(:),D2(:),GeoD(:,:)
Integer(4),allocatable::IE2(:),ITyp2(:),IDC(:)
Character(255),allocatable:: FilInp(:)

Integer(4) ITN2(MaxStr),iTypeStr(MaxStr,MaxStr),iir(1)

Character(30) buf30
Character(1) symb,symb3(4)

Integer(4) lRes(50),iRes(50)
Real(8) xRes(50)
Character(255) sRes(50)

Character(255) StrLabel,Str1,ToUpperCase
Real(8) E1value
Integer(4) iEword

Character(6) anm
Character(255) FilPat,FilGen
Character(132) line
Real(8) FreqX(3,MaxAt),Vibr1X(3,MaxAt)
Integer(4),allocatable::iStrings(:,:)
Logical lexist
Integer(4),allocatable::iEG1(:),iEG2(:),iEG3(:)
Real(8),allocatable::xEG1(:)


! Program for distinguishing of the cluster structures optimized with GA/ABC algorithms


! Get file names
FilDat='StrDiscr2.inp'
nrg=Command_Argument_Count()
If (nrg>0) Call GET_COMMAND_ARGUMENT(1,FilDat)
FilOut=FilDat
Call FileExtension(FilOut,'.out')

! Default thresholds
epsE=0.1d0
epsStr=0.03d0
epsPMOI=0.01d0
epsFP=1.d-4
Freq1min=-999999.d0

! Get job options
iDiscrMethodBeg=0
iDiscrMethodEnd=3
iE1value=0
iEtyp=2
iFilPat=0
iPrintUnit=6
ioutopen=0
iEinput=1
iShowPerPage=-1

If (ioutopen==0) Then
    Open(6,File=FilOut)
    Write(6,'(//''  *****  StrDiscr2  v.1.0 *****'')')
Endif

    iShowAll=0
    Nrow=1  ! Number of images in a row
    iUp2Down=1 ! 1/-1
    HorizontalShift=15.d0
    VerticalShift=15.d0

Open(3,File=FilDat)
Read(3,'(a255)')Str
ii=iGetOption('FILOUT',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    FilOut=sRes(1)
    Close(6, Status='DELETE')
    Open(6,File=FilOut)
    Write(6,'(//''  *****  StrDiscr2  v.1.0 *****'')')
    ioutopen=1
    iPrintUnit=6
Endif
ii=iGetOption('METHOD=STR',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    iDiscrMethodBeg=1
    iDiscrMethodEnd=1
Endif
ii=iGetOption('METHOD=PMOI',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    iDiscrMethodBeg=2
    iDiscrMethodEnd=2
Endif
ii=iGetOption('METHOD=FP',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    iDiscrMethodBeg=3
    iDiscrMethodEnd=3
Endif
ii=iGetOption('RBOND',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) RbondMax=xRes(1)
ii=iGetOption('EPSE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) epsE=xRes(1)
ii=iGetOption('EPSFP',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) epsFP=xRes(1)
ii=iGetOption('EPSSTR',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) epsStr=xRes(1)
ii=iGetOption('EPSPMOI',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) epsPMOI=xRes(1)
ii=iGetOption('STRLABEL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) StrLabel=sRes(1)
ii=iGetOption('EWORD',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEword=iRes(1)
ii=iGetOption('E1VALUE',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    E1value=xRes(1)
    iE1value=1
Endif
ii=iGetOption('EANALYSIS',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEtyp=iRes(1)
ii=iGetOption('NOGRAPH',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iGraphRecognition=0
ii=iGetOptionNoUcase('GEOTEMP',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    FilPat=sRes(1)
    iFilPat=1
Endif

ii=iGetOption('EINPUT',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEinput=iRes(1)
ii=iGetOption('EUNITS=AU',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEunits=1
ii=iGetOption('EUNITS=EV',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEunits=2
ii=iGetOption('EUNITS=KCAL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEunits=3
ii=iGetOption('EUNITS=KJ',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) iEunits=4

ii=iGetOption('SYMEPS',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) SymEps=xRes(1)
ii=iGetOption('FREQ1MIN',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Freq1min=xRes(1)
ii=iGetOption('SHOWALL',Str,iPrintUnit,nRes,lRes,iRes,xRes,sRes)
If (ii>0) Then
    iShowAll=1
    If (nRes>=1.and.iRes(1)>0) nRow=iRes(1)
    If (nRes>=2.and.iRes(2)>0) iUp2Down=iRes(2)
    If (nRes>=3.and.iRes(3)>0) HorizontalShift=xRes(3)
    If (nRes>=4.and.iRes(4)>0) VerticalShift=xRes(4)
    If (nRes>=5.and.iRes(5)>0) iShowPerPage=iRes(5)
Endif

Select Case(iEunits) 
Case(2) 
    Eunits=1.d0/27.21d0                  ! E in eV
Case(3) 
    Eunits=1.d0/(27.21d0*23.06d0)        ! E in kcal/mol
Case(4) 
    Eunits=1.d0/(27.21d0*23.06d0*4.184)  ! E in kJ/mol
Case Default
    Eunits=1.d0                          ! E in au by default
EndSelect




! Read input file name(s)

nFilInp=0
Open(1,File='tmp.tmp')
Do While(.not.EOF(3))
    Read(3,'(a255)')Str
    If (Len_Trim(Str)==0) Exit
    Str=AdjustL(Str)
    If (Str(1:1)=='!') Cycle
    nFilInp=nFilInp+1
    Write(1,'(a<Len_Trim(Str)>)')Trim(Str)
Enddo
Close(3)

Allocate(FilInp(nFilInp))
Rewind(1)
Do i=1,nFilInp
    Read(1,'(a255)')FilInp(i)
Enddo
Close(1,Status='DELETE')


Write(6,'(//'' Classification thresholds:'')')
Write(6,'('' Energy      '',f10.5)')epsE
Write(6,'('' Structure   '',f10.5)')epsStr
Write(6,'('' PMOI        '',f10.5)')epsPMOI
Write(6,'('' Fingerprints'',f10.5)')epsFP
Write(6,'('' SymGroup eps'',f10.5,'' A'')')symeps
If (rBondMax>0.d0) Then
    Write(6,'('' RbondMax    '',f10.5)')rBondMax
Else
    Write(6,'('' RbondMax not set -- will be used by default as 1.15*(Ra+Rb)'')')
Endif

Write(6,'(//'' Structures will be read in from files:'')')
Do i=1,nFilInp
    Write(6,'(i4,2x,a255)')i,FilInp(i)
Enddo
Write(6,'(/'' Structures will be selected from file using the label (case-insensitive): '',a<Len_Trim(StrLabel)>)')Trim(StrLabel)
Write(6,*)

! Read structures

istr=0
istr0=0
Etot=0.d0
au2kcal=27.21d0*23.06d0
isetup=0
AverDiam=0.d0
DiamMax=0.d0
DiamMin=999999.d0
idiamin=0
idiamax=0
Open(16,File='Structures-graphs.g6')
Open(11,File='Structures-vibdata.frq')
g6str=repeat(' ',Len(g6str))
Do inf=1,nFilInp
    Open(5,File=FilInp(inf))
    istrfil=0
    Do While (.not.EOF(5))
        If (INDEX(FilInp(inf),'.log')>0.or.INDEX(FilInp(inf),'.out')>0) Then       ! Read last geo, freq, vibr1 from G09/G16 log-file 
            FreqX=0.d0
            Vibr1X=0.d0
            Call ReadFreqs(5,istr,iopt,ierr,ierr2,ngeofound,netotfound,nfreq,nneg,EtotX,C,FreqX,Vibr1X)
            Str=Trim(FilInp(inf))
            Write(buf30,'(f30.8)')EtotX
            Str=Trim(Str)//'  '//Trim(AdjustL(buf30))
            If (iopt==2) Then
                Str=Trim(Str)//'  *opt  '
            Else
                Str=Trim(Str)//'  *end  '
            Endif
            If (nfreq>0) Str=Trim(Str)//'  *freq '
            If (nfreq>0.and.nneg==0) Then
                Str=Trim(Str)//'  *lm    '
            ElseIf (nfreq>0.and.nneg>0) Then
                Write(buf30,'(i30)')nneg
                buf30=AdjustL(buf30)
                Str=Trim(Str)//'  *nneg='//Trim(buf30)
            Endif
            ll=Len_Trim(Str)
            Write(Str1,'(<60-ll>x,3f10.2)')FreqX(1:3,1)
            Str=Trim(Str)//' '//Trim(Str1)
            If (Index(ToUpperCase(Str),Trim(StrLabel))==0) Cycle
            istr=istr+1
            StrLabelX(istr)=Str
            istr0=istr0+1
            istrfil=istrfil+1
            Freq(1:3,1:Numat,istr)=FreqX(1:3,1:Numat)
            Vibr1(1:3,1:Numat,istr)=Vibr1X(1:3,1:Numat)
            nvibr=nfreq
        Else                                        ! Read structures from xyz files with labels
            Read(5,'(a255)')Str
            Str=AdjustL(Str)
            If (Str(1:1)=='!') Cycle
            Str1=ToUpperCase(Str)
            If (INDEX(Str1,Trim(StrLabel))>0) Then
                istr=istr+1
                StrLabelX(istr)=Str
                istr0=istr0+1
                istrfil=istrfil+1
                Call SubString(Str,MaxSubStr,nSubStr,SubStr)
                buf30=SubStr(iEword)
                Read(buf30,'(f30.10)')EtotX
                Numat=0
                nfreq=0
                nvibr=0
                nneg=0
                line=repeat(' ',132)
                Do While (.not.EOF(5))
                    Read(5,'(a255)')Str
                    If (Len_Trim(Str)==0) Exit
                    isNaN1=0
                    isNaN2=0
                    isNaN3=0
                    Call SubString(Str,MaxSubStr,nSubStr,SubStr)
                    Do k=2,4
                        If (INDEX(SubStr(k),'NaN')>0) isNaN1=k
                        If (INDEX(SubStr(k),'***')>0) isNaN1=k
                    Enddo
                    Do k=5,nSubStr
                        If (INDEX(SubStr(k),'NaN')>0) isNaN2=k
                        If (INDEX(SubStr(k),'***')>0) isNaN2=k
                    Enddo                            
                    Do k=8,nSubStr
                        If (INDEX(SubStr(k),'NaN')>0) isNaN3=k
                        If (INDEX(SubStr(k),'***')>0) isNaN3=k
                    Enddo                            
                    If (isNaN1==0) Then
                        line=Trim(SubStr(1))//' '//Trim(SubStr(2))//' '//Trim(SubStr(3))//' '//Trim(SubStr(4))
                        Call ParseXYZ(Mode,NNuc,xx,yy,zz,line)
                        Numat=Numat+1
                        NA(Numat)=NNuc
                        C(1:3,Numat)=(/xx,yy,zz/)
                        If (nSubStr>=7.and.isNaN2==0) Then
                            Do k=1,3
                                nfreq=nfreq+1
                                buf30=SubStr(k+4)
                                Read(buf30,'(f30.10)')ff
                                Freq(k,Numat,istr)=ff
                                If (ff<0.d0) nneg=nneg+1
                            Enddo
                        Endif
                        If (nSubStr>=10.and.isNaN3==0) Then
                            Do k=1,3
                                nvibr=nvibr+1
                                buf30=SubStr(k+7)
                                Read(buf30,'(f30.10)')ff
                                Vibr1(k,Numat,istr)=ff
                            Enddo
                        Endif
                    Endif
                Enddo
            Else
                Cycle
            Endif
        Endif     
        EtotX=EtotX*Eunits      ! Convert E to au
        xNumat=Dble(Numat)
        If (iEinput==1) Then          !If E -> Etot
            Etot(1,istr)=EtotX
            Etot(2,istr)=(EtotX-E1value*xNumat)*au2kcal
            Etot(3,istr)=(EtotX/xNumat-E1value)*au2kcal
        ElseIf (iEinput==2) Then      ! If E--> Eb
            Etot(1,istr)=EtotX+E1value*Eunits*xNumat
            Etot(2,istr)=EtotX*au2kcal
            Etot(3,istr)=Etot(2,istr)/xNumat
        ElseIf (iEinput==3) Then      ! If E --> Eb/N
            Etot(2,istr)=EtotX*xNumat*au2kcal
            Etot(3,istr)=EtotX*au2kcal
            Etot(1,istr)=Etot(2,istr)/au2kcal+E1value*Eunits*xNumat
        Endif
        If (isNaN1>0) Then
            Write(6,'(i5,5x,3f20.8,'' Str discarded: NaN/coords'',10x,a<Len_Trim(StrLabelX(istr))>)')istr0,Etot(1:3,istr),Trim(StrLabelX(istr))
            istr=istr-1
            Cycle
        Endif
        If (isNaN2>0) Then
            Write(6,'(i5,5x,3f20.8,'' Str discarded: NaN/freqs '',10x,a<Len_Trim(StrLabelX(istr))>)')istr0,Etot(1:3,istr),Trim(StrLabelX(istr))
            istr=istr-1
            Cycle
        Endif
        If (isNaN3>0) Then
            Write(6,'(i5,5x,3f20.8,'' Str discarded: NaN/vibs  '',10x,a<Len_Trim(StrLabelX(istr))>)')istr0,Etot(1:3,istr),Trim(StrLabelX(istr))
            istr=istr-1
            Cycle
        Endif
        ff=Freq(1,1,istr)
        If (nfreq>0.and.ff<Freq1min) Then
            Write(6,'(i5,5x,3f20.8,'' Str discarded: Freq1:'',f10.2,4x,a<Len_Trim(StrLabelX(istr))>)')istr0,Etot(1:3,istr),ff,Trim(StrLabelX(istr))
            istr=istr-1
            Cycle
        Endif
        If (nneg>0) nfreq=-nfreq
        NFreqStr(istr)=nfreq
        NVibrStr(istr)=nvibr
        If (isetup==0) Then
            isetup=1
            Numat0=Numat
            Call SetupGraphs
            Write(6,'(/''    i Fil FilStr Str         Energy(1)           Energy(2)           Energy(3)  StrDiam  Connect       Graph    Graph No.  SymGrp      StrLabel in file'')')
        Endif
        If (Numat/=Numat0) Then
            Write(6,'('' ERROR! Number of atoms in the structure '',i5,'' is different from the previous structures'',i6)')istr
            Stop
        Endif
        Do i=1,Numat
            Amass(i)=AMS(NA(i))
        Enddo
        Call InertiaNew(2,1,1,Numat,C,Amass,CM,PMOI,PAxes,kRotTyp)
        Cstr(1:3,1:Numat,istr)=C(1:3,1:Numat)
        NAstr(1:Numat,istr)=NA(1:Numat)
        PM(1:3,istr)=PMOI(1:3)
        NAstr(1:Numat,istr)=NA(1:Numat)
        NumatStr(istr)=Numat
        Call GraphConnectivity(iconn)
        If (iGraphRecognition==1) Then
            Call DetermineGraphType(istr,igtype)
            IsoGraph(istr)=igType 
            If (igtype==-1) Call GraphString
        Else
            IsoGraph(istr)=-1
        Endif
        If (iConn>1) Then
            Write(6,'(i5,5x,3f20.8,f8.2'' Str discarded! Conn:'',i5,10x,a<Len_Trim(StrLabelX(istr))>)')istr0,Etot(1:3,istr),StrDiameter,iConn,Trim(StrLabelX(istr))
            istr=istr-1
            Cycle
        Endif
        Call SymGroup(kRotTyp,iSym,iSymGrp,NRotSym,SymGrp)
        SymGrpStr(istr)=SymGrp
        Write(6,'(i5,2i4,i5,3f20.8,f8.2,i5,12x,a10,i8,2x,a4,8x,a<Len_Trim(StrLabelX(istr))>)')istr0,inf,istrfil,istr,Etot(1:3,istr),StrDiameter,iConn,g6str(1:10),igtype,SymGrp,Trim(StrLabelX(istr))
        Write(16,'(a<Len_Trim(g6str)>)')Trim(g6str)
        AverDiam=AverDiam+StrDiameter
        If (StrDiameter>DiamMax) Then
            DiamMax=StrDiameter
            idiamax=istr
        Endif
        If (StrDiameter<DiamMin) Then
            DiamMin=StrDiameter
            idiamin=istr
        Endif
    Enddo

    Close(5)

Enddo
Close(16)

nstr=istr
AverDiam=AverDiam/Dble(nstr)
Write(6,'(/'' Structures found in file(s)  :'',i5)')istr0
Write(6,'( '' Structures to be analyzed    :'',i5)')istr
Write(6,'( '' Number of atoms in each str  :'',i5)')Numat0
Write(6,'( '' Average structure diameter, A:'',f10.2,'' (among'',i5,'' accepted structures)'')')AverDiam,istr
Write(6,'( '' Minimum structure diameter, A:'',f10.2,'' (among'',i5,'' accepted structures) istr:'',i8)')DiamMin,istr,idiamin
Write(6,'( '' Maximum structure diameter, A:'',f10.2,'' (among'',i5,'' accepted structures) istr:'',i8)')DiamMax,istr,idiamax



! Graph assignment
!If (Numat==2) FilCan='graph2c_canon.g6'
!If (Numat==3) FilCan='graph3c_canon.g6'
!If (Numat==4) FilCan='graph4c_canon.g6'
!If (Numat==5) FilCan='graph5c_canon.g6'
!If (Numat==6) FilCan='graph6c_canon.g6'
!If (Numat==7) FilCan='graph7c_canon.g6'
!If (Numat==8) FilCan='graph8c_canon.g6'
!If (Numat==9) FilCan='graph9c_canon.g6'
iCanonic=0
FilCan='Structures-canonic.g6'
Inquire(File=FilCan,EXIST=lexist)
If (lexist) Then
    iCanonic=1
    Write(6,'(//'' Classification on the basis of canonical graphs (canonization by McKay NAUTY method using the external program labelg)'')')
    Write(6,'(  '' Reading canonic g6 strings from file '',a<Len_Trim(FilCan)>,'' ...''\)')Trim(FilCan)
    Open(16,File=FilCan)
    lmax=0
    mstr=0
    Do While (.not.EOF(16))
        Read(16,'(a100)')g6str
        ll=Len_Trim(g6str)
        If (ll==0) Cycle
        If (ll>lmax) lmax=ll
        mstr=mstr+1
    Enddo
    Write(6,'(i10,'' g6 strings have been read in.'',a<Len_Trim(FilCan)>)')mstr,FilCan
    Allocate(iStrings(0:lmax,nstr))
    Rewind(16)
    is=0
    Do While (.not.EOF(16))
        Read(16,'(a<lmax>)')g6str(1:lmax)
        If (is==nstr) Exit
        is=is+1
        ll=Len_Trim(g6str)
        iStrings(0,is)=ll
        Do i=1,ll
            iStrings(i,is)=ICHAR(g6str(i:i))
        Enddo
    Enddo
    Close(16)    
    ns=is
    ic=0
    Do i=1,ns
        li=iStrings(0,i)
        If (li<0) Cycle
        ic=ic+1
        Write(6,'(/''Canonic graph type'',i4,2x,<li>a1,''  str:'',i4\)')ic,(CHAR(iStrings(jj,i)),jj=1,li),i
        IsoGraph(i)=ic
A1:      Do j=1,ns
            If (j==i) Cycle
            lj=iStrings(0,j)
            If (lj<0) Cycle A1
            If (li/=lj) Cycle A1
            Do k=1,lj
                If (iStrings(k,i)/=iStrings(k,j)) Cycle A1
            Enddo
            Write(6,'(i4\)')j
            iStrings(0,j)=-iStrings(0,j)
            IsoGraph(j)=ic
        Enddo A1
    Enddo
    Write(6,'(//''Classes (canonic graph types) found:'',i5)')ic
Else
    Write(6,'(//'' WARNING!  Important notice:'')')
    Write(6,'(  '' Can not perform classification on the basis of canonical graphs -- file with canonic g6 strings not found: '',a<Len_Trim(FilCan)>)')Trim(FilCan)
    Write(6,'(  '' To perform this analysis, process the generated file Structures-graph.g6 with the program labelg (./labelg Structures-graph.g6 Structures-canonic.g6)'')')
    Write(6,'(  '' and place its output file Structures-canonic.g6 into the current directory.'')')
Endif

    
!
! Calculation sorted energies
!
Allocate(E1(nstr))
Allocate(IE1(nstr),ITT1(nstr))

E1(1:nstr)=Etot(iEtyp,1:nstr)
Do i=1,nstr
    IE1(i)=i
Enddo

Call HeapSort1(nstr,E1,IE1)


!
! Calculate sorted interatomic distances
!
ndc=Numat*(Numat-1)/2
Allocate(DC(ndc,nstr),XDC(ndc))
Allocate(IE2(nstr),ITyp2(nstr),IDC(ndc))
Do i=1,nstr
    IE2(i)=i
Enddo
Do istr=1,nstr
    ij=0
    Do i=2,Numat
        Do j=1,i-1
            C(1:3,1:Numat)=Cstr(1:3,1:Numat,istr)
            rij=Distance(i,j,Numat,C)
            ij=ij+1
            XDC(ij)=rij
            IDC(ij)=ij
        Enddo
    Enddo
    Call HeapSort1(ndc,XDC,IDC)
    DC(1:ndc,istr)=XDC(1:ndc)
Enddo

!
! Main cycle
!
Do iDiscrMethod=iDiscrMethodBeg,iDiscrMethodEnd
Write(6,'(//80(''*''))')
Write(6,'(  80(''*''))')
Write(6,'(  80(''*''))')
ntypes=1
ITyp2=1

If (iDiscrMethod==0) Then   ! Discrimination by sorted interatomic distances
Write(6,'(/'' Discrimination by Energy('',i0,'').  Eps ='',1pe12.2,'' units of energy'')')iEtyp,epsE

eps=epsE
eps2=epsE*0.01d0
ittype=1
ITT1(1:nstr)=1
Do i=2,nstr
    If ((DABS(E1(i)-E1(i-1)))>eps) ittype=ittype+1
    ITT1(i)=ittype
    iTyp2(i)=ittype
Enddo
ntyp=ITT1(nstr)
ntypes=ntyp

Write(6,'(/''    i   Str     Energy(1)      Energy(2)      Energy(3) AnalyzedEnergy   ClassNo.         StrLabel in file'')')
Do i=1,nstr
    Write(6,'(2i5,4f15.6,i10,10x,a20)')i,IE1(i),Etot(1:3,IE1(i)),E1(i),ITT1(i),StrLabelX(IE1(i))
Enddo
!Write(6,'(/'' Classes distinguished by Energy('',i0,''):'',i0,i6)')iEtyp,ntyp

ElseIf (iDiscrMethod==1) Then   ! Discrimination by sorted interatomic distances
Write(6,'(''Discrimination by interatomic distances.''\)')
eps2=epsStr !0.03d0
    A: Do i=2,nstr
        XDC(1:ndc)=DC(1:ndc,i)
        Do j=1,i-1
            dxmax=0.d0
            Do k=1,ndc
                tmp=DABS((DC(k,j)-XDC(k))/XDC(k))
                If (tmp>dxmax) dxmax=tmp
            Enddo
            If (dxmax<epsStr) Then
                ITyp2(i)=ITyp2(j)
                Cycle A
            Endif
        Enddo
        ntypes=ntypes+1
        ITyp2(i)=ntypes
    Enddo A
ElseIf (iDiscrMethod==2) Then ! Discrimination by sorted PMOI
    Write(6,'(''Discrimination by PMOI values.''\)')
    eps2=epsPMOI    !0.01d0
    B: Do i=2,nstr
        PMOI(1:3)=PM(1:3,i)
        Do j=1,i-1
            dxmax=0.d0
            Do k=1,3
                tmp=DABS((PM(k,j)-PMOI(k))/PMOI(k))
                If (tmp>dxmax) dxmax=tmp
            Enddo
            If (dxmax<eps2) Then
                ITyp2(i)=ITyp2(j)
                Cycle B
            Endif
        Enddo
        ntypes=ntypes+1
        ITyp2(i)=ntypes
    Enddo B
ElseIf (iDiscrMethod==3) Then ! Discrimination by fingerprints
    eps2=epsFP
    Write(6,'('' Discrimination by Goedecker2013 fingerprints.  eps='',f8.5)')eps2
!    Write(6,'('' Fingerprints of all structures:'')')
    n4=Numat*4
    Allocate(D1(n4),D2(n4),GeoD(n4,nstr))
    Do i=1,nstr
        C(1:3,1:Numat)=Cstr(1:3,1:Numat,i)
        NA(1:Numat)=NAstr(1:Numat,i)
        Call GeoMetrics(0,Numat,NA,C,Score,D1)
        GeoD(1:n4,i)=D1(1:n4)
 !       Write(6,'(i5)')i
 !       Write(6,'(10g12.5)')D1(1:n4)
    Enddo    
    CC: Do i=2,nstr
        D1(1:n4)=GeoD(1:n4,i)
        Do j=1,i-1
            D2(1:n4)=GeoD(1:n4,j)
            D2=D2-D1
            Score=dot_product(D2,D2)/Dble(n4)
!!            Write(6,'(2i4,g25.6)')i,j,Score
            If (Score<eps2) Then
                ITyp2(i)=ITyp2(j)
                Cycle CC
            Endif
        Enddo
        ntypes=ntypes+1
        ITyp2(i)=ntypes
    Enddo CC
!    Write(6,*)
Endif
Write(6,'(/'' Classes found:'',i5)')ntypes

!
! Assignment structures to each class
!
iii=0
Do i=1,ntypes
    ii=0
    Do istr=1,nstr
        If (ITyp2(istr)==i) Then
            ii=ii+1
            iii=iii+1
            ITN2(i)=ii
            iTypeStr(i,ii)=istr
        Endif
    Enddo
Enddo

!
! Energy of classes
!
!Write(6,'(//'' Energies of different classes, ntypes='',i5)')ntypes
!Write(6,'(  ''  Typ Nstr      Eaver(1)          RMSE             Eaver(2)        RMSE         Eaver(3)        RMSE'')')
Eaver=0.d0
Estdv=0.d0
Do i=1,ntypes
    nti=ITN2(i)
    Do j=1,nti
        jj=iTypeStr(i,j)
        Eaver(1:3,i)=Eaver(1:3,i)+Etot(1:3,jj)
    Enddo
    Eaver(1:3,i)=Eaver(1:3,i)/Dble(nti)    ! Mean energy
    Do j=1,nti
        jj=iTypeStr(i,j)
        Estdv(1:3,i)=Estdv(1:3,i)+(Etot(1:3,jj)-Eaver(1:3,i))**2
    Enddo
    Estdv(1:3,i)=DSQRT(Estdv(1:3,i)/Dble(nti))
    Do k=1,3
        symb3(k)=' '
        If (DABS(Estdv(k,i)/Eaver(k,i))>0.005d0) symb3(k)='?'
    Enddo
!    Write(6,'(2i5,2f16.8,1x,a1,3x,2f12.4,1x,a1,3x,2f12.4,1x,a1)')i,nti,(Eaver(j,i),Estdv(j,i),symb3(j),j=1,3)
Enddo

!
! Sorting classes by energy
!
E1(1:ntypes)=Eaver(2,1:ntypes)
Call HeapSort1(ntypes,E1,IE1)

!
! Print original (not-transformed geometries - geometeries as they were optimized
!
If (iDiscrMethod>0) Then
    Write(6,'(//''Structural parameters of all structures. Nstr:'',i5,''   Ntypes:'',i5)')nstr,ntypes
    Write(6,'(  ''Structure  Str Type  Energy(2)   PMOI 1           2           3     Rij,A'')')
Endif
Open(8,File='Structures-original.xyz')
iii=0
ndcmax=Min(ndc,50)
Do i=1,ntypes
    ii=0
    Do istr=1,nstr
        If (ITyp2(istr)==i) Then
            ii=ii+1
            iii=iii+1
            If (iDiscrMethod>0) Write(6,'( ''@geo'',i0.4,2i5,f12.6,3f12.3,50f6.2)')istr,iii,i,Etot(iEtyp,istr),PM(1:3,istr),DC(1:50,istr)
            Write(8,'(/''@geo'',i0.4,2i5,f12.6,3g12.3,50f6.2)')istr,iii,i,Etot(iEtyp,istr),PM(1:3,istr),DC(1:ndcmax,istr)
            Do ia=1,Numat
                Write(8,'(i2,3f15.6)')NA(ia),Cstr(1:3,ia,istr)
            Enddo
        Endif
    Enddo
Enddo
Close(8)

! Reorient structures inside the class for better coincidence
Do i=1,ntypes
    i2i=ITN2(i)
    C(1:3,1:Numat)=Cstr(1:3,1:Numat,iTypeStr(i,1))
    NA(1:Numat)=NAstr(1:Numat,iTypeStr(i,1))
    Do j=2,i2i
        jj=iTypeStr(i,j)
        C1(1:3,1:Numat)=Cstr(1:3,1:Numat,jj)
        NA1(1:Numat)=NAstr(1:Numat,jj)
        dreor(0)=StrDiscrep(Numat,NA,NA1,C,C1)
        Do k=1,3
            Call RotXYZ(k,1,0,180.d0,0.d0,Numat,C1,C2)
            dreor(k)=StrDiscrep(Numat,NA,NA1,C,C2)
        Enddo
        Do k=4,6
            C2(1:3,1:Numat)=C1(1:3,1:Numat)
            C2(k-3,1:Numat)=-C1(k-3,1:Numat)
            dreor(k)=StrDiscrep(Numat,NA,NA1,C,C2)
        Enddo
        C2(1:3,1:Numat)=-C1(1:3,1:Numat)
        dreor(7)=StrDiscrep(Numat,NA,NA1,C,C2)
        iir=MinLoc(dreor)
        k=iir(1)-1
        If (k>=1.and.k<=3) Then
            Call RotXYZ(k,1,0,180.d0,0.d0,Numat,C1,C2)
            Cstr(1:3,1:Numat,jj)=C2(1:3,1:Numat)
        Endif
        If (k>=4.and.k<=6) Then
            C2(1:3,1:Numat)=C1(1:3,1:Numat)
            C2(k-3,1:Numat)=-C1(k-3,1:Numat)
            Cstr(1:3,1:Numat,jj)=C2(1:3,1:Numat)
        Endif
        If (k==7) Then
            C2(1:3,1:Numat)=-C1(1:3,1:Numat)
            Cstr(1:3,1:Numat,jj)=C2(1:3,1:Numat)
        Endif
    Enddo
Enddo

!
! Printing transformed structures -- structures transformed to be most similar to each another inside a class
!
Open(7,File='Structures-rotated.xyz')
Open(8,File='Structures-classes.xyz')
Write(7,'(''    istr    i  TypNo.  Typ  E(2)        E(3)     ii        PMOI(1)   PMOI(2)   PMOI(3)   dRij-sorted'')')
Write(8,'(''!First structures in their classes. Sorted by energy. Classification on the basis of ''\)')
If (iDiscrMethod==1) Write(8,'('' energies.'')')
If (iDiscrMethod==2) Write(8,'('' PMOI.'')')
If (iDiscrMethod==3) Write(8,'('' fingerprints.   E1value='',f16.8)')E1value
Write(8,'(''!   istr  Class iTyp         E(1)        E(2)        E(3)              PMOI(1)   PMOI(2)   PMOI(3)   dRij-sorted'')')
iii=0
Do it=1,ntypes
    i=IE1(it)
    ii=0
    Do istr=1,nstr
        If (ITyp2(istr)==i) Then
            ii=ii+1
            iii=iii+1
!            ITN2(i)=ii
!            iTypeStr(i,ii)=istr
!            Write(6,'( ''@geo'',i0.4,2i5,2f12.6,3i5,3g10.3,25f6.2)')istr,iii,i,Etot(2,istr),Etot(3,istr),ITyp2(istr),ITT1(istr),ii,PM(1:3,istr),DC(1:25,istr)
            Write(7,'(/''@geo'',i0.4,3i5,5x,3f12.6,i5,3g10.3,25f6.2)')istr,iii,it,i,Etot(1,istr),Etot(2,istr),Etot(3,istr),ii,PM(1:3,istr),DC(1:25,istr)
            Do ia=1,Numat
                Call SetAname(NA(ia),ia,anm)
                Write(7,'(a6,3f15.6)')anm,Cstr(1:3,ia,istr)
            Enddo
            If (ii==1) Then
                nfreq=IABS(nFreqStr(istr))
                nvibr=nVibrStr(istr)
                nneg=0
                If (nFreqStr(istr)<0) nneg=1
                buf5='     '
                If (nfreq>0) Then
                    buf5='*lm  '
                    If (nneg>0) buf5='*nneg'
                Endif
                Write(8,'(/''@geo'',i0.4,2i5,5x,3f12.6,2x,a5,3x,3g10.3,25f6.2)')istr,it,i,Etot(1,istr),Etot(2,istr),Etot(3,istr),buf5,PM(1:3,istr),DC(1:25,istr)
                Do ia=1,Numat
                    Call SetAname(NA(ia),ia,anm)
                    Write(8,'(a6,3f15.6\)')anm,Cstr(1:3,ia,istr)
                    If (nfreq>0) Write(8,'(5x,3f12.3\)')Freq(1:3,ia,istr)
                    If (nVibr>0) Write(8,'(5x,3f12.3\)')Vibr1(1:3,ia,istr)
                    Write(8,*)
                Enddo
                If (iFilPat==1.and.iDiscrMethod==3) Then
                    Open(9,File=FilPat)
                    Write(FilGen,'(''geo-c'',i0.3,''-s'',i0.5,''.gjf'')')it,istr
                    Open(10,File=FilGen)
                    Do While (.not.EOF(9))
                        Read(9,'(a255)')Str
                        Str1=ToUpperCase(Str)
                        If (INDEX(Str1,'%GEO%')>0) Then
                            Do ia=1,Numat
                                Call SetAname(NA(ia),ia,anm)
                                Write(10,'(a6,3f15.6)')anm,Cstr(1:3,ia,istr)
                            Enddo
                        Else
                            ll=Len_Trim(Str)
                            Write(10,'(a<ll>)')Trim(Str)
                        Endif
                    Enddo
                    Close(10)
                    Close(9)
                Endif
            Endif
        Endif
    Enddo
Enddo
Close(7)
Close(8)


!
! Final table
!
If (iDiscrMethod==0) Write(6,'(//'' Structure classes, discriminated by energy. Nstr:'',i5,''  Ntypes:'',i5)')nstr,ntypes
If (iDiscrMethod==1) Write(6,'(//'' Structure classes, discriminated by sorted interatomic distances. Nstr:'',i5,''  Ntypes:'',i5)')nstr,ntypes
If (iDiscrMethod==2) Write(6,'(//'' Structure classes, discriminated by PMOI. Nstr:'',i5,''  Ntypes:'',i5)')nstr,ntypes
If (iDiscrMethod==3) Write(6,'(//'' Structure classes, discriminated by fingerprints. Nstr:'',i5,''  Ntypes:'',i5)')nstr,ntypes
If (iDiscrMethod==0) Then
    Write(6,'(  '' Threshold for discrimination: '',g10.2,'' units of energy'',i2)')eps,iEtyp
Else
    Write(6,'(  '' Threshold for discrimination: '',f10.5,''% of each compared value.'')')eps2*100.d0
Endif
Write(6,'(  '' Classes with large energy scatter (SD) are marked with ? sign.'')')
Write(6,'(  '' Classes are sorted by energy.  SymGroup eps:'',f10.4,'' A'')')symeps
Write(6,'( /''Class iTyp Nstr Nbeg Nend    Eaver(1)    SD        Eaver(2)    SD        Eaver(3)    SD    IsoGraph   SymGrp   Structures:'')')
iend=0
Do it=1,ntypes
    i=IE1(it)
    i2i=ITN2(i)
    ibeg=iend+1
    iend=ibeg+i2i-1
    symb3(1:4)=' '
    If (Estdv(1,i)>1.0d0) symb3(1)='?'
    If (Estdv(2,i)>1.0d0) symb3(2)='?'
    If (Estdv(3,i)>0.1d0) symb3(3)='?'
    Do istr=1,nstr
        If (ITyp2(istr)==i) Exit
    Enddo
    If (iCanonic==1) Then
        igtype=IsoGraph(istr)
        Do j=1,i2i
            ijstr=iTypeStr(i,j)
            If (IsoGraph(ijstr)/=igtype) symb3(4)='?'
        Enddo
    Else
        igtype=-1
    Endif
    Write(6,'(5i5,3(f12.4,f8.4,1x,a1),i8,a1,2x,a4,2x,<i2i>i4)')it,i,ITN2(i),ibeg,iend,(Eaver(j,i),Estdv(j,i),symb3(j),j=1,3),igtype,symb3(4),SymGrpStr(istr),iTypeStr(i,1:i2i)
    If (symb3(4)=='?') Then
        If (iCanonic==0) Cycle
        Allocate(iEG1(i2i))
        k=0
        Do j=1,i2i
            ijstr=iTypeStr(i,j)
            itmp=IsoGraph(ijstr)
            kk=0
            Do l=1,k
                If (itmp==iEG1(l)) Then
                    kk=1
                    Exit
                Endif
            Enddo
            If (kk==0) Then
                k=k+1
                iEG1(k)=itmp
            Endif
        Enddo
        Write(6,'('' WARNING! Class'',i4,'' correpsonds to'',i4,'' canonic graph types:'')')it,k
        Do l=1,k
            itmp=iEG1(l)
            Write(6,'('' Graph '',i4,'' ->''\)')itmp
            Do j=1,i2i
                ijstr=iTypeStr(i,j)
                jtmp=IsoGraph(ijstr)
                If (jtmp==itmp) Write(6,'(i4\)')ijstr
            Enddo
            Write(6,*)
        Enddo
        !Write(6,*)
        Deallocate(iEG1)
    Endif
Enddo

Write(6,'( /''Class 1Str  SymGrp  Structure label in file.   SymGroup eps:'',f10.4,'' A'')')symeps
Do it=1,ntypes
    i=IE1(it)
    Do istr=1,nstr
        If (ITyp2(istr)==i) Then
            Write(6,'(2i5,2x,a4,4x,a<Len_Trim(StrLabelX(istr))>)')it,istr,SymGrpStr(istr),Trim(StrLabelX(istr))
            Exit
        Endif
    Enddo
Enddo

! Print g6 string for classes
If (iCanonic==1) Then
    Write(6,'( /''Class 1Str  SymGrp  Graph No. Canonic graph g6 string'')')
Do it=1,ntypes
    i=IE1(it)
    i2i=ITN2(i)
    Do istr=1,nstr
        If (ITyp2(istr)==i) Then
            igtype=IsoGraph(istr)
            Do j=1,i2i
                ijstr=iTypeStr(i,j)
                If (IsoGraph(ijstr)/=igtype) symb3(4)='?'
            Enddo
            ll=iStrings(0,istr)
            Write(6,'(2i5,2x,a4,4x,i5,5x,<ll>a1)')it,istr,SymGrpStr(istr),igtype,(CHAR(iStrings(jj,istr)),jj=1,ll)
            Exit
        Endif
    Enddo
Enddo
Else
    Write(6,'( /'' Canonic graph g6 file (Structures-canonic.g6) is missing. Canonic graph assignment skipped.'')')
Endif

! Show all structures at once on a single drawing
If (iShowAll==1) Then
    Open(9,File='ShowAll.xyz')
    NumatAll=ntypes*Numat
    Write(9,'(i8)')NumatAll
    Write(9,'(''! ''i5'' classes shown at once. Columns:'',i3\)')ntypes,nRow
    If (iUp2Down==1) Then
        Write(9,'(3x,''Direction:  Up->Down  Left->Right.''\)')
    Else
        Write(9,'(3x,''Direction:  Down->Up  Left->Right.''\)')
    Endif
    Write(9,'(3x,''Horizontal spacing:'',f6.2,'' Vertical spacing:'',f6.2)')HorizontalShift,VerticalShift
    ishow=-1
    ipage=0
    Do it=1,ntypes
        i=IE1(it)
        Do istr=1,nstr
            If (ITyp2(istr)==i) Then
                ishow=ishow+1
                If (ishowperpage>0.and.ishow/ishowperpage*ishowperpage==ishow) Then
                    ipage=ipage+1
                    Write(9,'(/''@page'',i0.2)')ipage
                Endif
                irow=ishow/nrow
                icol=ishow-(irow)*nrow
                !write(6,*)ishow,irow,icol
                dir=-Dble(iUp2Down)
                dx=Dble(icol)*HorizontalShift
                dy=Dble(irow)*VerticalShift*dir
                dz=0.d0
                Do ia=1,Numat
                    NA(ia)=NAstr(ia,istr)
                    C(1,ia)=Cstr(1,ia,istr)+dx
                    C(2,ia)=Cstr(2,ia,istr)+dy
                    C(3,ia)=Cstr(3,ia,istr)+dz
                Enddo
                !C=C*0.529177d0
                Do ia=1,Numat
                    Call SetAname(NA(ia),ia,Aname10)
                    Write(9,'(a10,3f15.6)')Aname10,C(1:3,ia)
                Enddo
                Exit
            Endif
        Enddo
    Enddo
    Close(9)
Endif

enddo

Close(11)

End
!***********************************************************************
Function StrDiscrep(Numat,NA,NA1,C,C1)
Implicit Real(8) (A-H,O-Z)

Real(8) C(3,Numat),C1(3,Numat)
Integer(4) NA(Numat),NA1(Numat),IA(Numat)

IA=0
StrDiscrep=0.d0
Do i=1,Numat
    rijmin=999999.d0
    jmin=0
    Do j=1,Numat
        If (NA(i)/=NA1(j)) Cycle
        If (IA(j)>0) Cycle
        rij=DSQRT((C(1,i)-C1(1,j))**2+(C(2,i)-C1(2,j))**2+(C(3,i)-C1(3,j))**2)
        If (rij<rijmin) Then
            rijmin=rij
            jmin=j
        Endif
    Enddo
    IA(jmin)=i
    StrDiscrep=StrDiscrep+rijmin
Enddo

End
