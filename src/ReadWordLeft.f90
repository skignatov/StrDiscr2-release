Subroutine ReadWordLeft(String,Str2)

Character(len=*) String,Str2
Character(len=Len(String)) Str1

! Subroutine ReadWordLeft reads a word (sequence of non-blank characters)
! from the beginning of Str1 and puts it to Str2

Str1=AdjustL(String)
l1=Len(Str1)
l2=Len(Str2)

If (l1==0.or.Len_Trim(Str1)==0) Then
	Str2(1:l2)=' '
	Return
Endif	

k=0
Do i=1,l1
	If (Str1(i:i)==' ') Exit
	k=k+1
	If (k>l2) Exit
	Str2(k:k)=Str1(i:i)
Enddo

If (k<l2) Str2(k+1:l2)=' '

Return
End
