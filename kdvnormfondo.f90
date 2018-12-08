!KORTEWEG-DEVRIES


!L'equazione risolta è quella di korteweg-deVries al caso normalizzato su fondovariabile
!in questa simulazione il fondo è con una funzione di -arctangente(scalino)
!Calcolo massimo,velocità 



program fondovariabile
integer i,icicli,isalva,i0,i1,i2,imax
real*8,allocatable, dimension (:) :: U,Unew,F,x,H,Xnew,vel
real*8,allocatable, dimension (:) :: xx_max,uu_max,vv_max,aa_max
character*30 filename
character nom1, nom2, nom0
real*8 :: pigreco,alpha,beta,Umax,xmax
!-------------------------------------variabili
common /parametri/ deltax,h0


pigreco=4.0*atan(1.0)

!------------------------Apro il file da cui prendo i dati

open(7,file='Dati.txt')
read(7,*) N
read(7,*) xmin
read(7,*) xmax
read(7,*) deltat
read(7,*) Ncicli
read(7,*) Nsalva
read(7,*) h0
read(7,*) c
read(7,*) a

close(1)
!---------------------------Ho preso i dati da file


allocate(U(0:N-1),Unew(0:N-1),x(0:N-1),F(0:N-1),H(0:N-1),Xnew(1:Ncicli),vel(1:Ncicli))
 deltax=(xmax-xmin)/N
 
 write(*,*) 'Valore di N pari a :',N
 write(*,*) 'Deltax uguale a :',deltax
 
 
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++ Condizioni iniziali
 
 !**********************************************************
 !----------------------------Alpha e Beta devono essere << 1
 alpha=a/h0 
 beta=((h0/c)**2)/6
 !**********************************************************
 
 Umax=0
 
 write(*,*) 'Inizializzazione'
 
 do i=0, N-1
 x(i)=xmin+i*deltax
 !write(*,*),x(i)
 
 U(i)=(dcosh(x(i))**(-2)) !controllare formula Whitham
 !print*,U(i)
 H(i)=0.2*datan(2*(x(i)-2.5))!****************************************tipo di fondo
  
end do
 !---------------------------------------------------------------- Fine condizioni Iniziali 
 
 
 
 write(*,*) 'Variabili Iniazializzate!'
 
 !...Apertura file salvataggio posizione picco
 open(20,file='pos_val_max.txt')
 
 
 !INIZIO CICLO
 !**********************************************************************
 do icicli=1, Ncicli
 !***********************************************************************
 
 
 !------------------------------------------------------start salva
 do isalva=0, Nsalva-1
 
 
 !+++++++++++++++++++++++++++++++++++++++++++++ inizio calcolo effettivo
 
  call calcolaeffe(U,F,N,H,alpha,beta)
  
  
   do i=0, N-1
     
      Unew(i)=U(i)+0.5*deltat*F(i)
   
   end do
   
  call calcolaeffe(Unew,F,N,H,alpha,beta)
  
  
   do i=0, N-1
  
     U(i)=U(i)+deltat*F(i)   
   end do
 
   Umax=U(0)
   do i=1,N-1
      if (U(i).gt.Umax) then
       Umax=U(i)
       imax=i
     end if
  
  end do
  
   x_max = x(imax)+.25*deltax*(U(imax+1)-U(imax-1))/(U(imax)-.5*(U(imax+1)+U(imax-1)))
   tempo=deltat*(isalva+(icicli-1)*Nsalva)
   write(20,*) tempo, x_max, Umax
 
  !+++++++++++++++++++++++++++++++++++++++++++++fine calcolo effettivo
 
 
 end do 
 !-------------------------------------------------------------------------------end salva
 write(*,*)'Il valore massimo di U è...................', Umax
 
 
 write(*,*) 'Ciclo numero',icicli
 
 !definiamo il nome dei file su cui scrivere, si utilizza il seguente metodo
 
 i0=  icicli/100
 i1=  (icicli-i0*100)/10
 i2=  icicli-i0*100-i1*10
 
 nom0=char(i0+ichar ('0'))
 nom1=char(i1+ichar ('0'))
 nom2=char(i2+ichar ('0'))
 
 
         filename='kdvnormform'//nom0//nom1//nom2//'.txt'
 
 write(*,*)'Inizio scrittura su File', filename
 write(*,*) x(icicli),'/', U(N/2), F(N/2)
 
 open(5,file=filename,form='formatted')
 do i=0, N-1
 write(5,300) x(i),U(i)
 end do
 close(5)
 
 
  
 
 write(*,*) 'Fine scrittura su File', filename
 
 !******************************************************************************
 end do
 !*******************************************************************************
 !FINE CICLO

 300 format(2e16.7)


 
 close(20)
 
 !---------------------------------------calcolo la velocità dell'onde onde

 
 
 deallocate(U,Unew,x,F,H,Xnew)
 
 
!...Aggiunta calcolo velocità picco

open(20,file='pos_val_max.txt')
allocate(xx_max(1:Ncicli*Nsalva),uu_max(1:Ncicli*Nsalva),vv_max(1:Ncicli*Nsalva),aa_max(1:Ncicli*Nsalva))
do i=1,Ncicli*Nsalva
read(20,*) tempo, xx_max(i), uu_max(i)
end do
close(20)
vv_max(1)=(xx_max(2)-xx_max(1))/deltat

do i=2,Ncicli*Nsalva-1
vv_max(i)=(xx_max(i+1)-xx_max(i-1))/(2*deltat)
aa_max(i)=(xx_max(i+1)-2*xx_max(i)+xx_max(i-1))/(deltat**2)
end do

vv_max(Ncicli*Nsalva)=(xx_max(Ncicli*Nsalva)-xx_max(Ncicli*Nsalva-1))/deltat

open(20,file='pos_val_max_vel_acc.txt')
do i=1,Ncicli*Nsalva
write(20,*) i*deltat, xx_max(i), uu_max(i), vv_max(i),aa_max(i)
end do
close(20)

deallocate(uu_max,vv_max,xx_max,aa_max)



 
 
end program fondovariabile




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ SUBROUTINE

subroutine calcolaeffe (U,F,N,H,alpha,beta)
real*8, dimension(0:N-1) :: U,F,GG,H
integer :: i
real*8 :: beta,alpha

common /parametri/ deltax,h0

do i=0,N-1
GG(i) = U(i)*(1.0+0.75*U(i)*alpha)
end do

do i=2, N-2


F(i)=-(GG(i+1)-GG(i-1))/(2*deltax)&
                - beta*(U(i+2)-2*U(i+1)+2*U(i-1)-U(i+2))/(2*deltax**3)&
                 + 0.25*U(i)*(H(i+1)-H(i-1))/deltax
    
   
 end do
 
 H(N-1)=H(0)
 
 F(N-1)=-(GG(0)-GG(N-2))/(2*deltax)&
                  - beta*(U(1)-2*U(0)+2*U(N-2)-U(N-3))/(2*deltax**3)&
                           +0.25*U(N-1)*(H(0)-H(N-2))/deltax
                    
 F(N-2)=-(GG(N-1)-GG(N-3))/(2*deltax)&
                  -beta*(U(0)-2*U(N-1)+2*U(N-3)-U(N-4))/(2*deltax**3)&
                           +0.25*U(N-2)*(H(N-1)-H(N-3))/deltax
 
 
 F(0)=-(GG(1)-GG(N-1))/(2*deltax)&
            - beta*((U(2)-2*U(1)+2*U(N-1)-U(N-2))/(2*deltax**3))&
                           +0.25*U(0)*(H(1)-H(N-1))/deltax

 F(1)=-(GG(2)-GG(0))/(2*deltax)&
             - beta*(U(3)-2*U(2)+2*U(0)-U(N-1))/(2*deltax**3)&
                           +0.25*U(1)*(H(2)-H(0))/deltax
  
 
 

end subroutine calcolaeffe

