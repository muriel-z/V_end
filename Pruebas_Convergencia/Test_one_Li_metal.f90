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
   logical,allocatable                     ::metal(:),mascara(:,:,:)
   !real(np)                                ::d1,d2,d3,dist2,div,res
   real(np)                                ::res,mres
   real(np)                                ::boxmin(3),boxmax(3),h(3)
   integer                                 ::i,l
   real(np),dimension(:,:),allocatable     :: r       ! Posiciones
   !real(np),dimension(:), allocatable      ::hist
   !real(np)                                ::dx,x
   integer                                 ::ri,rj,rk
   integer,allocatable                     :: ceros(:,:)
   !integer                                 ::ibin,nbin
   !character(len=6)                        ::aux
   !character(len=15)                       ::filename
   integer                                 :: count_CG,count_Li, nceros
   !character(len=2)                        ::Li,CG

   !Este programa tiene que agarrar una configuracion inicial y calcular el potencial
   !Del programa de python

   !Dimensiones de la caja (final distinta a la inicial)
   boxmin(:)=[0._np,0._np,0._np]
   boxmax(:)=[300._np,300._np,300._np]

   !Variables enteras
   nvx=100
   nvy=100
   nvz=100
   N_iter=2000

   h(1)=(boxmax(1)-boxmin(1))/(nvx-1)
   h(2)=(boxmax(2)-boxmin(2))/(nvy-1)
   h(3)=(boxmax(3)-boxmin(3))/(nvz-1)

   !Leo el archivo de posición final a partir del cual vamos a calcular el potencial
   open(14,File='pos_final.xyz')
   read(14,*) pnum
   close(14)
   write(*,*) 'pnum = ',pnum

   ! XXX: Prueba con 1 particula
   pnum=1
   
   allocate(r(pnum,3))
   allocate(sym(pnum))
   allocate(metal(pnum))
   allocate(m(pnum))
   count_Li=0
   count_CG=0

   ! TODO dejarla crecer ajustar pnum
   allocate(ceros(8*pnum,3))

   ! XXX: Prueba con 1 particula
   ! open(14,File='pos_final.xyz')
   ! read(14,*)
   ! read(14,*)
   ! do i=1,pnum
   !    read(14,*) sym(i),r(i,:),m(i)
   !    if(sym(i)=='L') then
   !       count_Li=count_Li+1
   !       metal(i)=.false.
   !    elseif(sym(i)=='C') then
   !       count_CG=count_CG+1
   !       metal(i)=.true.
   !    endif
   ! enddo
   ! close(14)

   ! XXX: Prueba con 1 particula
   sym(1)='C'
   count_CG=count_CG+1
   metal(1)=.true.
   r(1,:)=[150._dp,150._dp,150._dp]


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
                           
   ! Matriz booleana
   allocate(mascara(nvx,nvy,nvz))
                  
   ! Establece los indices a volver cero
   ! y crea una mascara donde estan
   ! TODO todo lo que implique metal, posiciones
   ! o atomos, meterlo aca dentro. Y sacarlo
   ! de la interfaz (e.g. metal)
   call dendritas(metal,r,V0,nceros,ceros,mascara)

   do l=1,N_iter

      ! Set ceros
      do i=1,nceros
         V0(ceros(i,1),ceros(i,2),ceros(i,3))=0._np
      enddo

      call step_pbc(V0,nvx,nvy,nvz,V,res,mres,mascara)

      ! TODO: Residuo entre V y V0
      if(mod(l,4)==0) then
         ! res=sum((V0(1:nvx,1:nvy,1:nvz)-V(:,:,:))**2)
         write(*,*) 'N_iter=',l,"Res",res/(nvx*nvy*nvz),"MRes", mres
      endif

      ! Swap pointers
      Vtmp=>V0
      V0(0:,0:,0:)=>V(:,:,:)
      V(0:,0:,0:)=>Vtmp(:,:,:)
   enddo

   write(*,*) 'Ya terminó de calcular el potencial'
   call salida('Vp.dat',r,V,V0)

! 7 FORMAT (A3)
! 10 FORMAT (3(3X,A15))
! 11 FORMAT (3(2X,ES17.9))

contains

   subroutine dendritas(metal,r,V,nceros,ceros,mascara)
      implicit none

      logical,dimension(:),intent(in)        :: metal
      integer ,intent(out)                   :: ceros(:,:)
      logical ,intent(out)                   :: mascara(:,:,:)
      real(np),intent(in)                    :: r(:,:)
      real(np)                               :: ip(3),ia,ia2,d2
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(inout) :: V
      real(np),parameter                     :: a=3.2_np

      integer,intent(out)                    :: nceros
      integer                                :: n,pnum,iia
      integer                                :: i,j,k

      ! Valores inicales
      mascara(:,:,:)=.false.
      nceros=0

      !Condicion de metal - Se asigna voltaje nulo a la posciones de la malla donde haya estructura de dendritas
      ! !$OMP PARALLEL DO PRIVATE(N,RI,RJ,RK)
      pnum=size(r,1)
      do n=1,pnum ! sobre las particulas metalicas
         !write(*,*) 'recorro las particulas internas n=', n
         if(metal(n)) then

            !Nodo inferior mas cercano
            ip(:)=(r(n,:)-boxmin(:))/h(:)+1

            !El radio de corte en 
            !XXX: Esto asume que h(1)=h(2)=h(3)
            ia=a/h(1)
            if(2*ia<1._np) then
               print *, "WARNING: mejora la resolución de la malla"
            endif

            iia=int(ia)
            ia2=ia*ia
 
            !Nodo inferior mas cercano
            ri=int(ip(1))
            rj=int(ip(2))
            rk=int(ip(3))
                  
            do i= ri-iia, ri+iia+1
              do j= rj-iia, rj+iia+1
                do k= rk-iia, rk+iia+1

                  d2=(ip(1)-i)**2+(ip(2)-j)**2+(ip(3)-k)**2
                  if(d2>ia2) cycle

                  nceros=nceros+1
                  ceros(nceros,1)=i
                  ceros(nceros,2)=j
                  ceros(nceros,3)=k

                  mascara(i,j,k)=.true.

                  V(i,j,k)=0._np

                enddo
              enddo
            enddo

         endif
      enddo

      print *, "We have this polos", nceros

   endsubroutine dendritas

   subroutine step_pbc(V0,nvx,nvy,nvz,V,res,mres,mask)
      implicit none

      integer,intent(in) :: nvx,nvy,nvz
      logical,intent(in) ::mask(:,:,:)
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(in) :: V0
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(out) :: V
      real(np),intent(out) :: res,mres
      real(np)             :: aux

      integer :: i,j,k

      res=0._np
      mres=0._np
      open(123,file='res',status='replace')
      !$omp parallel do private(i,j,k,aux) reduction(+:res) reduction(max: mres)
      do k=1,nvz
         do j=1,nvy
            do i=1,nvx
                
               ! Skip calculo residuos y potencial en lugares que son cero
               if(mask(i,j,k)) then
                 V(i,j,k)=0._np
                 cycle
               endif
                     
               V(i,j,k)=(V0(i+1,j,k)+V0(i-1,j,k)+V0(i,j+1,k)+V0(i,j-1,k)+V0(i,j,k+1)+V0(i,j,k-1))/6._np
               aux=(V0(i,j,k)-V(i,j,k))**2
               mres=max(mres,aux)
               res=res+aux
               !if(abs(V(i,j,k)-V0(i,j,k))>1.e-5_dp) print *, "WARNING"
            enddo

            ! Plano yz en el borde x superior
            V(nvx+1,j,k)=V(1,j,k)
            ! Plano yz en el borde x inferior
            V(0,j,k)=V(nvx,j,k)
         enddo

         ! Plano xz en el borde y superior
         V(1:nvx,nvy+1,k)=V(1:nvx,1,k)
         ! Plano xz en el borde y inverior
         V(1:nvx,0,k)=V(1:nvx,nvy,k)

         ! Las aristas no son necesarias porque la aproximacion del
         ! laplaciano que tenemos no las usa.
         ! V(0,0,k)=V(nvx,nvy,k)
         ! V(nvx+1,nvy+1,k)=V(1,1,k)
         ! V(0,nvy+1,k)=V(nvx,1,k)
         ! V(nvx+1,0,k)=V(1,nvy,k)
         
      enddo
      !$omp end parallel do
      close(123)

   endsubroutine step_pbc
               
   subroutine salida_V(archivo,V)
      character(*),intent(in)  :: archivo
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(in) :: V
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
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(in) :: V0
      real(np),dimension(0:nvx+1,0:nvy+1,0:nvz+1),intent(in) :: V

      real(np) :: Vpp(size(r,1))
      real(np) :: r3(3)
      real(np) :: div
      integer :: l,pnum

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
