program quantum_pressure 
implicit none

    ! all the unit is atomic unit

    integer, parameter                          :: dp=kind(1.0d0)
    real,parameter                              :: pi=3.14159265358979
    real,parameter                              :: eq=1d-7
    real,parameter                              :: c0=137.0359996287515

    type atoms
      integer                                   :: neu
      real(dp)                                  :: neu_real
      real(dp),dimension(3)                     :: posi
    end type

    type charge_grid
      real(dp),dimension(:,:,:),allocatable     :: density
      integer,dimension(3)                      :: grid
      real(dp),dimension(3,3)                   :: spaces
      real(dp)                                  :: volumegrid
      integer                                   :: natom
      real(dp),dimension(3)                     :: origin
      type(atoms),dimension(  :),allocatable    :: atom
    end type

    type(charge_grid)                           :: rho_charge
    type(charge_grid)                           :: t_charge
    type(charge_grid)                           :: tw_charge
    type(charge_grid)                           :: p_charge
    type(charge_grid)                           :: j_charge
    type(charge_grid)                           :: current_charge_x, current_charge_y, current_charge_z
    type(charge_grid)                           :: a_charge_x, a_charge_y, a_charge_z
    type(charge_grid)                           :: a2_charge
    type(charge_grid)                           :: lap_r_charge
    type(charge_grid)                           :: reduced_lap_charge
    type(charge_grid)                           :: prox_charge

    character(len=256)                          :: rho_file,t_file,tw_file,lap_file,j_file,output_file,output_file2,arg_tmp,version
    character(len=256)                          :: c1_file, c2_file, c3_file
    real(dp),dimension(:,:,:),allocatable       :: tf


    real(dp)                                    :: lap_coeff=1d0/4d0, mag_field
    integer                                     :: narg, equ

    version = 'Oct_2014-1'
    write(6,*) 'VERSION = ' , version

    narg =iargc()

    call getarg(1,rho_file)                             ! 1. density file

    call getarg(2,t_file)                               ! 2. energy density file (-5/6*(tao-tao^tf)+|grad n|**2/(48*n)

    call getarg(3,c1_file)                              ! 3. paramagnetic current file x 

    call getarg(4,c2_file)                              ! 3. paramagnetic current file y 

    call getarg(5,c3_file)                              ! 3. paramagnetic current file z 

    call getarg(6,arg_tmp); read(arg_tmp,*)mag_field    ! 6. Magnetic field strength (a.u.)

    call getarg(7,output_file)                          ! 7. output filename

    output_file2=trim(output_file)//'A_2.cube'

    !! read rho charge
    call read_dgrid_cube(rho_file,rho_charge)

    !! read t charge
    call read_dgrid_cube(t_file,t_charge)

    !! read current density
    call read_dgrid_cube(c1_file,current_charge_x)
    call read_dgrid_cube(c2_file,current_charge_y)
    call read_dgrid_cube(c3_file,current_charge_z)

    !! calculate A (vector potential)
    call calc_vector_potential(mag_field,a_charge_x,a_charge_y,a_charge_z,rho_charge)

    !! estimate weak field term
    call calc_vector_square(a_charge_x,a_charge_y,a_charge_z,rho_charge,a2_charge)

    !! write cube file
    call write_cube(output_file2,a2_charge)

    !! calculate p^kq with normalization
    call calc_pkq_oct_14_v1(rho_charge,t_charge,a_charge_x,a_charge_y,a_charge_z,                   &
                            current_charge_x,current_charge_y,current_charge_z,p_charge)

    !! write cube file
    call write_cube(output_file,p_charge)

    !! finish
    write(6,*) 'Finished!'

    !
    !
    !
    !
    !
    !



contains 


    SUBROUTINE read_dgrid_cube(filename,charges)
        type(charge_grid)                                      :: charges
        character(len= * ),intent(in)                          :: filename
  
        integer                                                :: i, j, k
        character(len=256)                                     :: title
  
        open(unit=22,File=trim(filename),form='formatted')
  
        !- unit of bohr
        read(22,*) title
        read(22,*) title
        read(22,*) charges%natom, charges%origin(:)
  
        allocate( charges%atom(charges%natom) )
  
        read(22,*) charges%grid(1), charges%spaces(1,:)
        read(22,*) charges%grid(2), charges%spaces(2,:)
        read(22,*) charges%grid(3), charges%spaces(3,:)
        do i=1,charges%natom
            read(22,*) charges%atom(i)%neu, charges%atom(i)%neu_real, charges%atom(i)%posi(:)
        end do
  
        if (.not. allocated(charges%density))                                           &
            allocate(charges%density( charges%grid(1),charges%grid(2),charges%grid(3))  )
  
        read(22,*) (((charges%density(i,j,k),k=1,charges%grid(3)),j=1,charges%grid(2)),i=1,charges%grid(1))
  
        charges%volumegrid=charges%spaces(1,1)*charges%spaces(2,2)*charges%spaces(3,3)
  
        close(22)

    end SUBROUTINE




    SUBROUTINE calc_p_charge(rho_charge,t_charge,tw_charge,current_charge,p_charge,tf)
        type(charge_grid)                                    :: rho_charge
        type(charge_grid),intent(inout)                      :: t_charge
        type(charge_grid),intent(inout)                      :: tw_charge
        type(charge_grid)                                    :: current_charge
        type(charge_grid),intent(out)                        :: p_charge
        real(dp),dimension(:,:,:),allocatable                :: tf

        
        if (check_dimension(rho_charge,t_charge,tw_charge) ) then
            if (allocated(p_charge%density)) call deallo_charge(p_charge)
            p_charge=copy_cube_head(rho_charge)
            if (allocated(tf)) deallocate(tf)
            allocate(              tf(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)))
            tf=tf_charge(rho_charge%density)
     !      call deallo_charge(rho_charge)
            allocate(p_charge%density(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)) )
        else
            write(6,*) 'Dimension wrong for read files!'
            stop
        end if

        p_charge%density=1d0/6d0*tw_charge%density-5d0/6d0*(t_charge%density-tf)

     !  if (allocated(t_charge%density)) call deallo_charge(t_charge)
     !  if (allocated(tw_charge%density)) call deallo_charge(tw_charge)

    end SUBROUTINE




    SUBROUTINE calc_p_charge_current(rho_charge,t_charge,tw_charge,current_charge,p_charge,tf)
        type(charge_grid)                                    :: rho_charge
        type(charge_grid),intent(inout)                      :: t_charge
        type(charge_grid),intent(inout)                      :: tw_charge
        type(charge_grid),intent(inout)                      :: current_charge
        type(charge_grid),intent(out)                        :: p_charge
        real(dp),dimension(:,:,:),allocatable                :: tf

        
        if (check_dimension(rho_charge,t_charge,tw_charge) ) then
            if (allocated(p_charge%density)) call deallo_charge(p_charge)
            p_charge=copy_cube_head(rho_charge)
            if (allocated(tf)) deallocate(tf)
            allocate(              tf(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)))
            tf=tf_charge(rho_charge%density)
     !      call deallo_charge(rho_charge)
            allocate(p_charge%density(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)) )
        else
            write(6,*) 'Dimension wrong for read files!'
            stop
        end if

        p_charge%density=1d0/6d0*tw_charge%density-5d0/6d0*(t_charge%density-tf)-current_charge%density/(3d0*rho_charge%density)

     !  if (allocated(t_charge%density)) call deallo_charge(t_charge)
     !  if (allocated(tw_charge%density)) call deallo_charge(tw_charge)

    end SUBROUTINE




    SUBROUTINE pk_charge(p_charge,tf,zero_proc,one)
        type(charge_grid),intent(inout)                      :: p_charge
        real(dp),dimension(:,:,:)                            :: tf
        real(dp),intent(in)                                  :: zero_proc
        real(dp),intent(in)                                  :: one

        real(dp),dimension(size(tf,1),size(tf,2),size(tf,3)) :: tf_tmp
        integer                                              :: i, j, k
        real(dp)                                             :: eq_tmp

        tf_tmp=0d0
        tf=tf*2d0/3d0
        eq_tmp=1d-6**(5d0/3d0)*3d0/10d0*(3d0*pi**2)**(2d0/3d0)
        forall(i=1:size(tf,1),j=1:size(tf,2),k=1:size(tf,3),abs(tf(i,j,k))< eq_tmp*2d0/3d0)     &
            tf_tmp(i,j,k)=one+zero_proc/sqrt(1d0+zero_proc**2)

        forall(i=1:size(tf,1),j=1:size(tf,2),k=1:size(tf,3),abs(tf(i,j,k))>=eq_tmp*2d0/3d0)     &
            tf_tmp(i,j,k)=one+p_charge%density(i,j,k)/tf(i,j,k)/sqrt(1d0+(p_charge%density(i,j,k)/tf(i,j,k))**2)

        ! should make this expression better
        p_charge%density=tf_tmp
 

    end SUBROUTINE





!   SUBROUTINE write_radial_2charge(p_charge,charge_1,charge_2,tf,zero_proc,one,uni,scale_factor1,scale_factor2)
!       type(charge_grid),intent(inout)                      :: p_charge
!       real(dp),dimension(:,:,:)                            :: charge_1, charge_2
!       real(dp),dimension(:, :, :)                          :: tf
!       real(dp),intent(in)                                  :: zero_proc
!       real(dp),intent(in)                                  :: one
!       real(dp),intent(in)                                  :: scale_factor1, scale_factor2
!       integer,intent(in)                                   :: uni

!       real(dp),dimension(:),allocatable                    :: radial_r
!       real(dp),dimension(size(tf,1),size(tf,2),size(tf,3)) :: distance, charge_3, tf_tmp
!       integer                                              :: i, j, k, length
!       integer,dimension(3)                                 :: center
!       real(dp)                                             :: eq_tmp

!       distance=0d0
!       tf_tmp=tf*2d0/3d0
!       charge_3=charge_1*scale_factor1+charge_2*scale_factor1
!       eq_tmp=eq**(5d0/3d0)*3d0/10d0*(3d0*pi**2)**(2d0/3d0)
!       center=find_center(p_charge%grid, p_charge%spaces)
!       length=minval( p_charge%grid-center )-2
!       if (allocated(radial_r)) deallocate(radial_r)
!       allocate( radial_r(length+1) )

!       forall(i=1:p_charge%grid(1),j=1:p_charge%grid(2),k=1:p_charge%grid(3))    &
!                     distance(i,j,k)= ((i-center(1))*p_charge%spaces(1,1))**2   +&
!                                      ((j-center(2))*p_charge%spaces(2,2))**2   +&
!                                      ((k-center(3))*p_charge%spaces(3,3))**2
!       ! along first vector, orthom box
!       j=1
!       do i=center(2),center(2)+length
!           if (abs(tf(center(1),i,center(3))) < eq_tmp*2d0/3d0) radial_r(j)=one+zero_proc/sqrt(1d0+zero_proc**2)
!           if (abs(tf(center(1),i,center(3))) >=eq_tmp*2d0/3d0) radial_r(j)=one+charge_3(center(1),i,center(3))/tf_tmp(center(1),i,center(3))/sqrt(1d0+(charge_3(center(1),i,center(3))/tf_tmp(center(1),i,center(3)))**2)
!           j=j+1
!       end do

!       write(uni,'(2E30.12)') (sqrt(distance(center(1),i,center(3))), radial_r(i-center(2)+1), i=center(2),center(2)+length)

!       ! along ab
!       j=1
!       do i=center(1),center(1)+length
!           if (abs(tf(i,i,center(3))) < eq_tmp*2d0/3d0) radial_r(j)=one+zero_proc/sqrt(1d0+zero_proc**2)
!           if (abs(tf(i,i,center(3))) >=eq_tmp*2d0/3d0) radial_r(j)=one+charge_3(i,i,center(3))/tf_tmp(i,i,center(3))/sqrt(1d0+(charge_3(i,i,center(3))/tf_tmp(i,i,center(3)))**2)
!           j=j+1
!       end do

!       write(uni+1,'(2E30.12)') (sqrt(distance(i,i,center(3))), radial_r(i-center(1)+1), i=center(1),center(1)+length)

!       ! along abc
!       j=1
!       do i=center(1),center(1)+length
!           if (abs(tf(i,i,i)) < eq_tmp*2d0/3d0) radial_r(j)=one+zero_proc/sqrt(1d0+zero_proc**2)
!           if (abs(tf(i,i,i)) >=eq_tmp*2d0/3d0) radial_r(j)=one+charge_3(i,i,i)/tf_tmp(i,i,i)/sqrt(1d0+(charge_3(i,i,i)/tf_tmp(i,i,i))**2)
!           j=j+1
!       end do

!       write(uni+2,'(2E30.12)') (sqrt(distance(i,i,i)), radial_r(i-center(1)+1), i=center(1),center(1)+length)


!   end SUBROUTINE




    SUBROUTINE write_bare_radial(p_charge,density,zero_proc,one,uni,scale_factor)
        type(charge_grid),intent(inout)                      :: p_charge
        real(dp),dimension(:, :, :)                          :: density
        real(dp),intent(in)                                  :: zero_proc
        real(dp),intent(in)                                  :: one
        real(dp),intent(in)                                  :: scale_factor
        integer,intent(in)                                   :: uni

        real(dp),dimension(:),allocatable                    :: radial_r
        real(dp),dimension(size(density,1),size(density,2),size(density,3)) :: distance
        integer                                              :: i, j, k, length
        integer,dimension(3)                                 :: center
        real(dp)                                             :: eq_tmp

        distance=0d0
        p_charge%density=p_charge%density*scale_factor
        eq_tmp=eq**(5d0/3d0)*3d0/10d0*(3d0*pi**2)**(2d0/3d0)
        center=find_center(p_charge%grid, p_charge%spaces)
        length=minval( p_charge%grid-center )-2
        if (allocated(radial_r)) deallocate(radial_r)
        allocate( radial_r(length+1) )

        forall(i=1:p_charge%grid(1),j=1:p_charge%grid(2),k=1:p_charge%grid(3))    &
                      distance(i,j,k)= ((i-center(1))*p_charge%spaces(1,1))**2   +&
                                       ((j-center(2))*p_charge%spaces(2,2))**2   +&
                                       ((k-center(3))*p_charge%spaces(3,3))**2
        ! along first vector, orthom box
        j=1
        do i=center(2),center(2)+length
            radial_r(j)=density(center(1),i,center(3))
            j=j+1
        end do

        write(uni,'(2E30.12)') (sqrt(distance(center(1),i,center(3))), radial_r(i-center(2)+1), i=center(2),center(2)+length)

        ! along ab
        j=1
        do i=center(1),center(1)+length
            radial_r(j)=density(i,i,center(3))
            j=j+1
        end do

        write(uni+1,'(2E30.12)') (sqrt(distance(i,i,center(3))), radial_r(i-center(1)+1), i=center(1),center(1)+length)

        ! along abc
        j=1
        do i=center(1),center(1)+length
            radial_r(j)=density(i,i,i)
            j=j+1
        end do

        write(uni+2,'(2E30.12)') (sqrt(distance(i,i,i)), radial_r(i-center(1)+1), i=center(1),center(1)+length)


    end SUBROUTINE





!   SUBROUTINE write_radial(p_charge,tf,zero_proc,one,uni,scale_factor)
!       type(charge_grid),intent(inout)                      :: p_charge
!       real(dp),dimension(:, :, :)                          :: tf
!       real(dp),intent(in)                                  :: zero_proc
!       real(dp),intent(in)                                  :: one
!       real(dp),intent(in)                                  :: scale_factor
!       integer,intent(in)                                   :: uni

!       real(dp),dimension(:),allocatable                    :: radial_r
!       real(dp),dimension(size(tf,1),size(tf,2),size(tf,3)) :: distance
!       integer                                              :: i, j, k, length
!       integer,dimension(3)                                 :: center
!       real(dp)                                             :: eq_tmp

!       distance=0d0
!       p_charge%density=p_charge%density*scale_factor
!       eq_tmp=eq**(5d0/3d0)*3d0/10d0*(3d0*pi**2)**(2d0/3d0)
!       center=find_center(p_charge%grid, p_charge%spaces)
!       length=minval( p_charge%grid-center )-2
!       if (allocated(radial_r)) deallocate(radial_r)
!       allocate( radial_r(length+1) )

!       forall(i=1:p_charge%grid(1),j=1:p_charge%grid(2),k=1:p_charge%grid(3))    &
!                     distance(i,j,k)= ((i-center(1))*p_charge%spaces(1,1))**2   +&
!                                      ((j-center(2))*p_charge%spaces(2,2))**2   +&
!                                      ((k-center(3))*p_charge%spaces(3,3))**2
!       ! along first vector, orthom box
!       j=1
!       do i=center(2),center(2)+length
!           if (abs(tf(center(1),i,center(3))) < eq_tmp*2d0/3d0) radial_r(j)=one+zero_proc/sqrt(1d0+zero_proc**2)
!           if (abs(tf(center(1),i,center(3))) >=eq_tmp*2d0/3d0) radial_r(j)=one+p_charge%density(center(1),i,center(3))/(2d0/3d0*tf(center(1),i,center(3)))/sqrt(1d0+(p_charge%density(center(1),i,center(3))/(2d0/3d0*tf(center(1),i,center(3))))**2)
!           j=j+1
!       end do

!       write(uni,'(2E30.12)') (sqrt(distance(center(1),i,center(3))), radial_r(i-center(2)+1), i=center(2),center(2)+length)

!       ! along ab
!       j=1
!       do i=center(1),center(1)+length
!           if (abs(tf(i,i,center(3))) < eq_tmp*2d0/3d0) radial_r(j)=one+zero_proc/sqrt(1d0+zero_proc**2)
!           if (abs(tf(i,i,center(3))) >=eq_tmp*2d0/3d0) radial_r(j)=one+p_charge%density(i,i,center(3))/(2d0/3d0*tf(i,i,center(3)))/sqrt(1d0+(p_charge%density(i,i,center(3))/(2d0/3d0*tf(i,i,center(3))))**2)
!           j=j+1
!       end do

!       write(uni+1,'(2E30.12)') (sqrt(distance(i,i,center(3))), radial_r(i-center(1)+1), i=center(1),center(1)+length)

!       ! along abc
!       j=1
!       do i=center(1),center(1)+length
!           if (abs(tf(i,i,i)) < eq_tmp*2d0/3d0) radial_r(j)=one+zero_proc/sqrt(1d0+zero_proc**2)
!           if (abs(tf(i,i,i)) >=eq_tmp*2d0/3d0) radial_r(j)=one+p_charge%density(i,i,i)/(2d0/3d0*tf(i,i,i))/sqrt(1d0+(p_charge%density(i,i,i)/(2d0/3d0*tf(i,i,i)))**2)
!           j=j+1
!       end do

!       write(uni+2,'(2E30.12)') (sqrt(distance(i,i,i)), radial_r(i-center(1)+1), i=center(1),center(1)+length)


!   end SUBROUTINE





    FUNCTION find_center(num,spaces)

      integer,dimension(3)                      :: find_center, num
      real(dp),dimension(3:3)                   :: spaces

      find_center=(/num(1)/2+1, num(2)/2+1, num(3)/2+1 /)

    end FUNCTION





    SUBROUTINE comp_j2_rho(rho_charge, j_charge, p_charge)
        type(charge_grid),intent(inout)                      :: rho_charge
        type(charge_grid),intent(inout)                      :: j_charge
        type(charge_grid),intent(inout)                      :: p_charge
       

        if (allocated(p_charge%density)) call deallo_charge(p_charge)
        p_charge=copy_cube_head(rho_charge)
        allocate(p_charge%density(size(rho_charge%density,1),     &
                                  size(rho_charge%density,2),     &
                                  size(rho_charge%density,3)) )

        p_charge%density=j_charge%density/rho_charge%density

    end SUBROUTINE





    SUBROUTINE pq_charge(p_charge,lap_charge,t_charge,rho_charge,tf,lap_coeff)
        type(charge_grid),intent(inout)                      :: lap_charge
        type(charge_grid),intent(inout)                      :: t_charge
        type(charge_grid),intent(inout)                      :: rho_charge
        type(charge_grid),intent(inout)                      :: p_charge
        real(dp),           optional                         :: lap_coeff
        real(dp),dimension(:,:,:),allocatable                :: tf

        if (check_dimension(rho_charge,t_charge,lap_charge) ) then
            if (allocated(p_charge%density)) call deallo_charge(p_charge)
            p_charge=copy_cube_head(rho_charge)
            if (allocated(tf)) deallocate(tf)
            allocate(              tf(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)))
            tf=tf_charge(rho_charge%density)
     !      call deallo_charge(rho_charge)
            allocate(p_charge%density(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)) )
        else
            write(6,*) 'Dimension wrong for read files!'
            stop
        end if

        
        if (.not. present(lap_coeff)) lap_coeff=1d0/4d0
        p_charge%density=2d0/3d0*(t_charge%density-tf)-lap_coeff*lap_charge%density

     !  if (allocated(t_charge%density)) call deallo_charge(t_charge)
     !  if (allocated(lap_charge%density)) call deallo_charge(lap_charge)

    end SUBROUTINE





    SUBROUTINE pq_charge_current(p_charge,lap_charge,t_charge,rho_charge,j_charge,tf,lap_coeff)
        type(charge_grid),intent(inout)                      :: lap_charge
        type(charge_grid),intent(inout)                      :: t_charge
        type(charge_grid),intent(inout)                      :: rho_charge
        type(charge_grid),intent(inout)                      :: p_charge
        type(charge_grid),intent(inout)                      :: j_charge
        real(dp),           optional                         :: lap_coeff
        real(dp),dimension(:,:,:),allocatable                :: tf

        if (check_dimension(rho_charge,t_charge,lap_charge) ) then
            if (allocated(p_charge%density)) call deallo_charge(p_charge)
            p_charge=copy_cube_head(rho_charge)
            if (allocated(tf)) deallocate(tf)
            allocate(              tf(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)))
            tf=tf_charge(rho_charge%density)
     !      call deallo_charge(rho_charge)
            allocate(p_charge%density(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)) )
        else
            write(6,*) 'Dimension wrong for read files!'
            stop
        end if

        
        ! unit of current from nA/T to a.u. (1 a.u. = 28.17909 nA/T)
        if (.not. present(lap_coeff)) lap_coeff=1d0/4d0
        p_charge%density=2d0/3d0*(t_charge%density-tf)-lap_coeff*lap_charge%density-j_charge%density /(3d0*rho_charge%density)

     !  if (allocated(t_charge%density)) call deallo_charge(t_charge)
     !  if (allocated(lap_charge%density)) call deallo_charge(lap_charge)

    end SUBROUTINE





    SUBROUTINE pq_charge_exempt_tf(p_charge,lap_charge,t_charge,rho_charge,tf,lap_coeff)
        type(charge_grid),intent(inout)                      :: lap_charge
        type(charge_grid),intent(inout)                      :: t_charge
        type(charge_grid),intent(inout)                      :: rho_charge
        type(charge_grid),intent(inout)                      :: p_charge
        real(dp),           optional                         :: lap_coeff
        real(dp),dimension(:,:,:),allocatable                :: tf

        if (check_dimension(rho_charge,t_charge,lap_charge) ) then
            if (allocated(p_charge%density)) call deallo_charge(p_charge)
            p_charge=copy_cube_head(rho_charge)
            if (allocated(tf)) deallocate(tf)
            allocate(              tf(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)))
            tf=tf_charge(rho_charge%density)
     !      call deallo_charge(rho_charge)
            allocate(p_charge%density(size(rho_charge%density,1),     &
                                      size(rho_charge%density,2),     &
                                      size(rho_charge%density,3)) )
        else
            write(6,*) 'Dimension wrong for read files!'
            stop
        end if

        
        if (.not. present(lap_coeff)) lap_coeff=1d0/4d0
        p_charge%density=2d0/3d0*(t_charge%density)-lap_coeff*lap_charge%density

     !  if (allocated(t_charge%density)) call deallo_charge(t_charge)
     !  if (allocated(lap_charge%density)) call deallo_charge(lap_charge)

    end SUBROUTINE




    SUBROUTINE calc_pkq_oct_14_v1(rho_charge,t_charge,a_x,a_y,a_z,j_x,j_y,j_z,p_charge)
        type(charge_grid),intent(inout)                      :: rho_charge
        type(charge_grid),intent(inout)                      :: t_charge
        type(charge_grid),intent(inout)                      :: a_x, a_y, a_z
        type(charge_grid),intent(inout)                      :: j_x, j_y, j_z
        type(charge_grid),intent(inout)                      :: p_charge

        integer                                              :: i,j,k

        if (allocated(p_charge%density)) call deallo_charge(p_charge)
        p_charge=copy_cube_head(rho_charge)
        allocate(p_charge%density(size(rho_charge%density,1),     &
                                  size(rho_charge%density,2),     &
                                  size(rho_charge%density,3)) )

        p_charge%density=0d0

        forall(i=1:rho_charge%grid(1),j=1:rho_charge%grid(2),k=1:rho_charge%grid(3))  
            p_charge%density(i,j,k)=t_charge%density(i,j,k)                                          &
                + ( a_x%density(i,j,k)*j_x%density(i,j,k)                                            &
                +   a_y%density(i,j,k)*j_y%density(i,j,k)                                            &
                +   a_z%density(i,j,k)*j_z%density(i,j,k) )*2.d0/3d0/c0                              &
                + ( a_x%density(i,j,k)*a_x%density(i,j,k)                                            &
                +   a_x%density(i,j,k)*a_x%density(i,j,k)                                            &
                +   a_x%density(i,j,k)*a_x%density(i,j,k) )*rho_charge%density(i,j,k)/3.d0/c0**2      
        end forall

        call d_by_tf(p_charge,rho_charge,1d0,-1d0)

    end SUBROUTINE




    ! this function will modify the value of in_charge
    SUBROUTINE d_by_tf(in_charge,rho_charge,one,zero_proc)
        type(charge_grid),intent(inout)                      :: in_charge
        type(charge_grid),intent(in   )                      :: rho_charge
        real(dp),intent(in)                                  :: one
        real(dp),intent(in)                                  :: zero_proc

        real(dp),dimension(:,:,:),allocatable                :: tf
        real(dp)                                             :: eq_tmp
        integer                                              :: i,j,k



        if (allocated(tf)) deallocate(tf)
        allocate(              tf(size(rho_charge%density,1),     &
                                  size(rho_charge%density,2),     &
                                  size(rho_charge%density,3)))
        tf=tf_charge(rho_charge%density)

        tf=tf*2d0/3d0
        eq_tmp=1d-6**(5d0/3d0)*3d0/10d0*(3d0*pi**2)**(2d0/3d0)
        forall(i=1:size(tf,1),j=1:size(tf,2),k=1:size(tf,3),abs(tf(i,j,k))< eq_tmp*2d0/3d0)     &
            in_charge%density(i,j,k)=0.5d0*(one+zero_proc/sqrt(1d0+zero_proc**2))

        forall(i=1:size(tf,1),j=1:size(tf,2),k=1:size(tf,3),abs(tf(i,j,k))>=eq_tmp*2d0/3d0)     &
            in_charge%density(i,j,k)=0.5d0*(one+in_charge%density(i,j,k)                        &
                                    /tf(i,j,k)/sqrt(1d0+(in_charge%density(i,j,k)/tf(i,j,k))**2))

        if (allocated(tf)) deallocate(tf)


    end SUBROUTINE





    SUBROUTINE calc_vector_square(a_x,a_y,a_z,rho_charge,a2_charge)
        type(charge_grid),intent(in   )                      :: rho_charge
        type(charge_grid),intent(inout)                      :: a_x,a_y,a_z
        type(charge_grid),intent(inout)                      :: a2_charge

        integer                                              :: i,j,k


        if (allocated(a2_charge%density)) call deallo_charge(a2_charge)
        a2_charge=copy_cube_head(rho_charge)
        allocate(a2_charge%density(size(rho_charge%density,1),     &
                                   size(rho_charge%density,2),     &
                                   size(rho_charge%density,3)) )


        forall(i=1:rho_charge%grid(1),j=1:rho_charge%grid(2),k=1:rho_charge%grid(3))             &
            a2_charge%density(i,j,k)=rho_charge%density(i,j,k)                                   &
                * ( a_x%density(i,j,k)*a_x%density(i,j,k)                                        &
                +   a_y%density(i,j,k)*a_y%density(i,j,k)                                        &
                +   a_z%density(i,j,k)*a_z%density(i,j,k) ) / (3d0*c0**2)

!       call d_by_tf(a2_charge,rho_charge,1d0,-1d0)

    end SUBROUTINE





    SUBROUTINE calc_vector_potential(mag_field_z,a_x,a_y,a_z,rho_charge)
        type(charge_grid),intent(inout)                      :: a_x,a_y,a_z
        type(charge_grid),intent(in   )                      :: rho_charge
        real(dp),intent(in)                                  :: mag_field_z

        integer                                              :: num_x, num_y, num_z,i,j,k
        real(dp)                                             :: sx, sy, sz
        real(dp),dimension(3)                                :: coord

        if (allocated(a_x%density)) call deallo_charge(a_x)
        if (allocated(a_y%density)) call deallo_charge(a_y)
        if (allocated(a_z%density)) call deallo_charge(a_z)
        allocate(a_x%density(size(rho_charge%density,1),     &
                             size(rho_charge%density,2),     &
                             size(rho_charge%density,3)) )
        allocate(a_y%density(size(rho_charge%density,1),     &
                             size(rho_charge%density,2),     &
                             size(rho_charge%density,3)) )
        allocate(a_z%density(size(rho_charge%density,1),     &
                             size(rho_charge%density,2),     &
                             size(rho_charge%density,3)) )

        num_x=size(rho_charge%density,1)
        num_y=size(rho_charge%density,2)
        num_z=size(rho_charge%density,3)
        sx=rho_charge%spaces(1,1)
        sy=rho_charge%spaces(2,2)
        sz=rho_charge%spaces(3,3)

        write(6,*) 'ok'
        forall(i=1:num_x,j=1:num_y,k=1:num_z)
            a_x%density(i,j,k)=(j*sy+rho_charge%origin(2))*mag_field_z*(-1d0)*0.5d0
            a_y%density(i,j,k)=(i*sx+rho_charge%origin(1))*mag_field_z*( 1d0)*0.5d0
            a_z%density(i,j,k)=(k*sz+rho_charge%origin(3))*0d0
        end forall

        write(6,*) 'ik'

    end SUBROUTINE





    SUBROUTINE comp_approx(prox_charge,rho_charge,tw_charge,lap_charge)
        type(charge_grid),intent(inout)                      :: lap_charge
        type(charge_grid),intent(inout)                      :: tw_charge
        type(charge_grid),intent(inout)                      :: rho_charge
        type(charge_grid),intent(inout)                      :: prox_charge

        if (allocated(prox_charge%density)) call deallo_charge(prox_charge)
        prox_charge=copy_cube_head(rho_charge)
        if (allocated(tf)) deallocate(tf)
        allocate(              tf(size(rho_charge%density,1),     &
                                  size(rho_charge%density,2),     &
                                  size(rho_charge%density,3)))
        tf=tf_charge(rho_charge%density)
     !  call deallo_charge(rho_charge)
        allocate(prox_charge%density(size(rho_charge%density,1),     &
                                     size(rho_charge%density,2),     &
                                     size(rho_charge%density,3)) )
        prox_charge%density=tf+1d0/9d0*tw_charge%density+1d0/6d0*lap_charge%density

     !  if (allocated(tw_charge%density)) call deallo_charge(tw_charge)
     !  if (allocated(lap_charge%density))call deallo_charge(lap_charge)

    end SUBROUTINE




    SUBROUTINE comp_diff(p_charge,t_charge,prox_charge)
        type(charge_grid)                                    :: p_charge
        type(charge_grid),intent(in   )                      :: t_charge
        type(charge_grid),intent(in   )                      :: prox_charge
        
        if (allocated(p_charge%density)) call deallo_charge(p_charge)
        p_charge=copy_cube_head(t_charge)
        allocate(p_charge%density(size(rho_charge%density,1),     &
                                  size(rho_charge%density,2),     &
                                  size(rho_charge%density,3)) )

        p_charge%density=t_charge%density-prox_charge%density

    end SUBROUTINE





    SUBROUTINE comp_reduced_lap(r_lap_charge,lap_charge,tf)
        type(charge_grid),intent(inout)                      :: lap_charge
        type(charge_grid),intent(inout)                      :: r_lap_charge
        real(dp),dimension(:,:,:)                            :: tf

        if (allocated(r_lap_charge%density)) call deallo_charge(r_lap_charge)
        r_lap_charge=copy_cube_head(lap_charge)
        allocate(r_lap_charge%density(size(lap_charge%density,1),     &
                                      size(lap_charge%density,2),     &
                                      size(lap_charge%density,3)) )

        r_lap_charge%density=lap_charge%density/(40d0/3d0*tf)

    end SUBROUTINE





    SUBROUTINE deallo_charge(charge)
        type(charge_grid)                                    :: charge

        if (allocated(charge%atom)) deallocate(charge%atom)
        if (allocated(charge%density)) deallocate(charge%density)

    end SUBROUTINE




    FUNCTION tf_charge(density)
        real(dp),dimension(:,:,:)                            :: density
        real(dp),dimension(size(density,1),size(density,2),size(density,3))  :: tf_charge

        tf_charge=3d0/10d0*(3d0*pi**2)**(2d0/3d0)*density**(5d0/3d0)

    end FUNCTION




    FUNCTION check_dimension(charge1,charge2,charge3)
        type(charge_grid),intent(in)                         :: charge1, charge2, charge3
        logical                                              :: check_dimension
        check_dimension=.true.
    end FUNCTION




    FUNCTION copy_cube_head(charge)
        type(charge_grid),intent(in)                         :: charge
        type(charge_grid)                                    :: copy_cube_head
      
        integer                                              :: i

        copy_cube_head%natom=charge%natom
        copy_cube_head%volumegrid=charge%volumegrid
        copy_cube_head%spaces=charge%spaces
        copy_cube_head%grid=charge%grid
        copy_cube_head%origin=charge%origin

        allocate( copy_cube_head%atom(charge%natom) )

        do i=1,charge%natom
            copy_cube_head%atom(i)%neu=charge%atom(i)%neu
            copy_cube_head%atom(i)%neu_real=charge%atom(i)%neu_real
            copy_cube_head%atom(i)%posi=charge%atom(i)%posi
        end do

    end FUNCTION




    SUBROUTINE write_cube(filename,charge)
        character(len=*),intent(in)                          :: filename
        type(charge_grid),intent(inout)                      :: charge

        integer                                              :: i, j, k

        open(21,file=filename,form='formatted',access='append')

        if (allocated(charge%atom)) then
            write(21,*) 'Quantum pressure space grid,', ' version =', trim(version)
            write(21,*) 'Cube format'
    
            write(21,'(I8,3F12.6)') charge%natom, charge%origin
            write(21,'(I8,3F12.6)') charge%grid(1),charge%spaces(1,:)
            write(21,'(I8,3F12.6)') charge%grid(2),charge%spaces(2,:)
            write(21,'(I8,3F12.6)') charge%grid(3),charge%spaces(3,:)
            do i=1,charge%natom
                write(21,'(I8,4F12.6)') charge%atom(i)%neu, charge%atom(i)%neu_real, charge%atom(i)%posi(:)
            end do
        else
            write(6,*) 'Atom info not found'
        end if

        if (allocated(charge%density)) then
            write(21,'(6E18.8)') (((charge%density(i,j,k),k=1,charge%grid(3)),j=1,charge%grid(2)),i=1,charge%grid(1))
        else
            write(6,*) 'density is empty!'
        end if

        close(21)
    end SUBROUTINE




end program


!!!!
