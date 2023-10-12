program V_test_end

   use precision

   implicit none

   integer                                 ::nvx,nvy,nvz !numero de grilla
   integer                                 ::pnum,N_iter !numero de particulas, numero de iteraciones
   real(np),dimension(:,:,:,:),allocatable,target :: Vstorage
   real(np),dimension(:,:,:),pointer       ::V,V0,Vtmp        !voltaje nuevo y anterior
   !real(np),dimension(:,:,:),allocatable   ::Vx,Vy,Vz,V0x,V0y,V0z !componentes
   !real(np)                                ::dv
   real(np)                                ::Vtop !Vtop es el valor del voltaje en la tapa
   real(np),dimension(:),allocatable       ::m   !Voltaje asociado a cada particula,masa
   character,dimension(:),allocatable      ::sym         !indica si la particula es litio libre o congelado
   logical,dimension(:),allocatable        ::metal
   !real(np)                                ::d1,d2,d3,dist2,div,res
   real(np)                                ::res
   real(np)                                ::boxmin(3),boxmax(3),h(3)
   integer                                 ::i,l
   real(np),dimension(:,:),allocatable     :: r       ! Posiciones
   !real(np),dimension(:), allocatable      ::hist
   !real(np)                                ::dx,x
   integer                                 ::ri,rj,rk
   !integer                                 ::ibin,nbin
   !character(len=6)                        ::aux
   !character(len=15)                       ::filename
   integer                                 :: count_CG,count_Li
   !character(len=2)                        ::Li,CG

   !Este programa tiene que agarrar una configuracion inicial y calcular el potencial
   !Del programa de python

   !Dimensiones de la caja (final distinta a la inicial)
   boxmin(:)=[0._np,0._np,0._np]
   boxmax(:)=[300._np,300._np,300._np]

   !Variables enteras
   nvx=50
   nvy=50
   nvz=50
   N_iter=2000

   h(1)=(boxmax(1)-boxmin(1))/(nvx-1)
   h(2)=(boxmax(2)-boxmin(2))/(nvy-1)
   h(3)=(boxmax(3)-boxmin(3))/(nvz-1)

   !Leo el archivo de posición final a partir del cual vamos a calcular el potencial
   open(14,File='pos_final.xyz')
   read(14,*) pnum
   close(14)
   write(*,*) 'pnum = ',pnum

   allocate(r(pnum,3))
   allocate(sym(pnum))
   allocate(metal(pnum))
   allocate(m(pnum))
   count_Li=0
   count_CG=0


   open(14,File='pos_final.xyz')
   read(14,*)
   read(14,*)
   do i=1,pnum
        read(14,*) sym(i),r(i,:),m(i)
        if(sym(i)=='L') then
            count_Li=count_Li+1
            metal(i)=.false.
        elseif(sym(i)=='C') then
        count_CG=count_CG+1
        metal(i)=.true.
        endif
   enddo
   close(14)


   write(*,*) 'count_Li = ',count_Li

   write(*,*) 'count_CG = ',count_CG

   write(*,*) 'Ya leyó las posiciones finales'

   !Variables del potencial
   allocate(Vstorage(0:nvx+1,0:nvy+1,0:nvz+1,2))

   V0(0:,0:,0:)=>Vstorage(:,:,:,1)
   V(0:,0:,0:)=>Vstorage(:,:,:,2)

   write(*,*) 'size V',size(V,dim=1)
   write(*,*) 'size V0',size(V0,dim=3)

   !Condicion inicial - construyo el gradiente de voltaje en la dirección z
   Vtop=10._np
   do i=0,nvz+1
      V0(0:,0:,i)=real(i,dp)/(nvz+1)*Vtop
   enddo

   V(0:,0:,0:)=V0(0:,0:,0:) ! TODO: refine
   ! V(1:nvx,1:nvy,1:nvz)=V0(2:nvx+1,2:nvy+1,2:nvz+1)

   call salida('Vini.dat',r,V,V0)

   !Doy el valor de la matriz del voltaje en t=0.Resuelta por Laplace

   !Itero para encontrar la solucion suave

   res=sum((V0(0:,0:,0:)-V(0:,0:,0:))**2)
   write(*,*) 'Initial residue: ',res

   do l=1,N_iter
      call dendritas(metal,r,V0)

      call step_pbc(V0,nvx,nvy,nvz,V,res)

      ! TODO: Residuo entre V y V0
      if(mod(l,2)==0) then
         ! res=sum((V0(1:nvx,1:nvy,1:nvz)-V(:,:,:))**2)
         write(*,*) 'N_iter=',l,"Res",res
      endif

      ! Swap pointers
      Vtmp=>V0
      V0(0:,0:,0:)=>V(:,:,:)
      V(0:,0:,0:)=>Vtmp(:,:,:)
   enddo

   write(*,*) 'Ya terminó de calcular el potencial'

   print *, "RESS1SS",(V0(24,24,24)-V(24,24,24))**2
   call salida('Vp.dat',r,V,V0)

! 7 FORMAT (A3)
! 10 FORMAT (3(3X,A15))
! 11 FORMAT (3(2X,ES17.9))

contains

   subroutine dendritas(metal,r,V)
      implicit none

      logical,dimension(:),intent(in)        :: metal
      real(np),intent(in)                    :: r(:,:)
      real(np),dimension(:,:,:),intent(out)  :: V

      integer :: n,pnum


      !Condicion de metal - Se asigna voltaje nulo a la posciones de la malla donde haya estructura de dendritas
      !$OMP PARALLEL DO PRIVATE(N,RI,RJ,RK)
      do n=1,pnum ! sobre las particulas metalicas
         !write(*,*) 'recorro las particulas internas n=', n
         if(metal(n)) then
            !write(*,*)'n=',n
            ri=int((r(n,1)-boxmin(1))/h(1))+1
            !write(*,*) 'ri=', ri
            rj=int((r(n,2)-boxmin(2))/h(2))+1
            !write(*,*) 'rj =', rj
            rk=int((r(n,3)-boxmin(3))/h(3))+1
            !write(*,*) 'rk=', rk

            V(ri,rj,rk)=0._np
         endif
      enddo
   endsubroutine dendritas

   subroutine step(V0,nvx,nvy,nvz,V,res)
      implicit none

      integer,intent(in) :: nvx,nvy,nvz
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(in) :: V0
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(out) :: V
      real(np),intent(out) :: res

      integer :: i,j,k

      res=0._np
      !$omp parallel do private(i,j,k) reduction(+:res)
      do k=1,nvz
         do j=1,nvy
            do i=1,nvx
               V(i,j,k)=(V0(i+1,j,k)+V0(i-1,j,k)+V0(i,j+1,k)+V0(i,j-1,k)+V0(i,j,k+1)+V0(i,j,k-1))/6._np
               res=res+((V0(i,j,k)-V(i,j,k))**2)
               !if(abs(V(i,j,k)-V0(i,j,k))>1.e-5_dp) print *, "WARNING"
            enddo
         enddo
      enddo
      !$omp end parallel do
   endsubroutine step

   subroutine step_pbc(V0,nvx,nvy,nvz,V,res)
      implicit none

      integer,intent(in) :: nvx,nvy,nvz
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(in) :: V0
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(out) :: V
      real(np),intent(out) :: res

      integer :: i,j,k

      res=0._np
      open(123,file='res',status='replace')
      !$omp parallel do private(i,j,k) reduction(+:res)
      do k=1,nvz
         do j=1,nvy
            do i=1,nvx
               V(i,j,k)=(V0(i+1,j,k)+V0(i-1,j,k)+V0(i,j+1,k)+V0(i,j-1,k)+V0(i,j,k+1)+V0(i,j,k-1))/6._np
               !Esto esta mal pero es masomenos lo que tengo que hcaer para que no me haga este calculo
               if((V0(i,j,k)-V(i,j,k))**2.ge.9.8070400018248165_np) then
                 write(123,*) i,j,k,(V0(i,j,k)-V(i,j,k))**2
                 cycle
               endif
               res=res+((V0(i,j,k)-V(i,j,k))**2)
               !if(abs(V(i,j,k)-V0(i,j,k))>1.e-5_dp) print *, "WARNING"
            enddo
            V(nvx+1,j,k)=V(1,j,k)
            V(0,j,k)=V(nvx,j,k)
         enddo
         V(1:nvx,nvy+1,k)=V(1:nvx,1,k)
         V(1:nvx,0,k)=V(1:nvx,nvy,k)

         ! Cruzo aristas
         V(0,0,k)=V(nvx,nvy,k)
         V(nvx+1,nvy+1,k)=V(1,1,k)

         ! Cruzo aristas
         V(0,nvy+1,k)=V(nvx,1,k)
         V(nvx+1,0,k)=V(1,nvy,k)

      enddo
      !$omp end parallel do
      close(123)
   endsubroutine step_pbc

   subroutine salida_V(archivo,V)
      character(*),intent(in)  :: archivo
      real(np),dimension(:,:,:),intent(in) :: V
      integer :: i,j,k

      open(216,file=archivo,status='replace')
      print *, 'entro'
      do k=1,nvz
         do j=1,nvy
            do i=1,nvx
               write(216,'(3(i0,x),e15.7)') i,j,k,V(i,j,k)
            enddo
            write(216,*)
         enddo
         write(216,*)
         write(216,*)
      enddo
      print *, 'salio'
      close(216)

   endsubroutine

   subroutine salida(archivo,r,V,V0)
      character(*),intent(in)  :: archivo
      real(np),intent(in)      :: r(:,:)
      real(np),dimension(:,:,:),intent(in) :: V
      real(np),dimension(:,:,:),intent(in) :: V0

      real(np) :: Vpp(size(r,1))
      real(np) :: r3(3)
      real(np) :: div
      integer :: l,pnum

   print *, "RESS2SS",(V0(25,25,25)-V(25,25,25))**2
      call salida_V('V.dat',V)
      call salida_V('V0.dat',V0)

      pnum=size(r,1)

      open(16,file=archivo,status='replace')
      write(16,8) 'l','r(l,1)','r(l,2)','r(l,3)','Vp'

      Vpp(:)=0._np

      !Este loop deberia ser solamente sobre las particulas Li, o sea las que estan libres en el electrolito
      do l=1,pnum
         div=0._np
         ! write(*,*) real(l,np)/real(pnum,np)*100._np

         r3(:)=r(l,:)
         Vpp(l)=trilinear_interpolation(r3,V,boxmin,h)

         write(16,9) real(l,np),r(l,1),r(l,2),r(l,3),Vpp(l)
      enddo
      close(16)

8     FORMAT(5(3X,A15))
9     FORMAT(5(2X,ES17.9))

   endsubroutine

   function trilinear_interpolation(r,f,boxmin,h) result(interp_val)
      real(np),intent(in) :: r(3),f(:,:,:)
      real(np),intent(in) :: boxmin(3),h(3)
      real(np) :: interp_val

      integer :: i,j,k,i1,j1,k1
      real(np) :: x,y,z,dx,dy,dz

      i=int((r(1)-boxmin(1))/h(1))+1
      j=int((r(2)-boxmin(2))/h(2))+1
      k=int((r(3)-boxmin(3))/h(3))+1

      ! Esto deberia devolver error al salirse de la grilla
      ! if (i < 1) then
      !   i = 1
      ! elseif (i > size(f, 1) - 1) then
      !   i = size(f, 1) - 1
      ! endif
      !
      ! if (j < 1) then
      !   j = 1
      ! elseif (j > size(f, 2) - 1) then
      !   j = size(f, 2) - 1
      ! endif
      !
      ! if (k < 1) then
      !   k = 1
      ! elseif (k > size(f, 3) - 1) then
      !   k = size(f, 3) - 1
      ! endif

      i1=i+1
      j1=j+1
      k1=k+1

      x=(r(1)-boxmin(1))/h(1)-float(i-1)
      y=(r(2)-boxmin(2))/h(2)-float(j-1)
      z=(r(3)-boxmin(3))/h(3)-float(k-1)

      dx=1.0-x
      dy=1.0-y
      dz=1.0-z

      interp_val=dx*dy*dz*f(i,j,k)+x*dy*dz*f(i1,j,k) &
                  +dx*y*dz*f(i,j1,k)+x*y*dz*f(i1,j1,k) &
                  +dx*dy*z*f(i,j,k1)+x*dy*z*f(i1,j,k1) &
                  +dx*y*z*f(i,j1,k1)+x*y*z*f(i1,j1,k1)

   endfunction trilinear_interpolation

endprogram V_test_end
