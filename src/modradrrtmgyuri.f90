module modradrrtmgyuri
use modraddata

implicit none

contains
!__________________________________________________________
!__________________________________________________________
!__________________________________________________________
!__________________________________________________________
!Important Routine
  ! subroutine findGLQPoints(n_GLQ_clear, GLQ_points_clear, GLQ_weights_clear, n_clear, &
    ! n_GLQ_cloudtop, GLQ_points_cloudtop, GLQ_weights_cloudtop, n_clouds, &
    ! n_classes, n_in_class, class_size, total_amount_GLQ_points, GLQ_index_all, &
    ! original_clear_LWP_indexes, original_cloudtop_LWP_indexes)
  subroutine findGLQPoints()
  
    use modglobal, only: imax, jmax, kmax, i1, j1, k1, kind_rb, grav, zf, zh
    use modfields, only: ql0, qt0
    use modradrrtmgyuriroutines
    
    implicit none
    
    integer :: i,j,k,inverse_k, n, nbin
    integer :: im, jm, ksounding !Changed This (PIER_QV)
    
    !acceleration
    integer :: n_RT != nx * ny / n_RT_Ratio    !Actual number of Radtiative transfer calculations
    
    !Gauss-Legendre Quadrature
    integer :: N_g                              !Loops through GLQ_points
    integer :: x_index                            !Index of GLQ x_i
    integer :: temp_i, temp_j                        !Temporary placeholders for i and j in sorting and ordering
    integer :: GLQ_counter                          !counter necessary for full list of GLQ indexes
    integer,allocatable,dimension(:,:):: GLQ_clear_LWP_indexes      !The indexes necessary for the GLQ of the cloudless LWP
    integer,allocatable,dimension(:,:,:):: GLQ_cloudtop_LWP_indexes       !The indexes necessary for the GLQ of the cloudtop LWP, extra axis for the classes
    real(kind=kind_rb) :: binwidth, GLQ_val, valuewidth
    integer, dimension(:) :: temploc(1)
    
    !!Contained in modraddata
    !----------------------------------------------
    ! integer :: n_GLQ_cloudtop, n_GLQ_clear                  !Amount of points for GLQ (cloudtop is per class)
    ! integer :: total_amount_GLQ_points                    !Total amount of GLQ points (n_GLQ_clear + n_GLQ_cloudtop*n_classes)
    ! integer,allocatable,dimension(:,:):: GLQ_index_all          !All GLQ indexes in a single array starting with cloudless and appending the first clouded class after being followed by second clouded etc.
    ! integer,allocatable,dimension(:,:):: original_clear_LWP_indexes       !original indexes of cloudless_LWP_ordered
    ! integer,allocatable,dimension(:,:,:):: original_cloudtop_LWP_indexes    !original indexes of the sorted LWP
    ! real(kind=kind_rb),allocatable,dimension(:) :: GLQ_points_clear, GLQ_weights_clear  !GLQ values cloudless
    ! real(kind=kind_rb),allocatable,dimension(:,:) :: GLQ_points_cloudtop, GLQ_weights_cloudtop  !GLQ values cloudtop, extra axis for the classes
    
    ! integer,allocatable,dimension(:) :: n_in_class                 !Array that contains the amount of clouds in a certain class
    ! integer :: class_size              !Amount of clouds in individual class
    ! integer :: n_clouds, n_clear          !number of collums with clouds and number of clear collumns
    ! integer :: n_classes_initial                !maximum number of cloudtop altitude classes
    ! real(kind=kind_rb)      :: cloud_threshold             !for the definition of a clouded collumn
    ! real(kind=kind_rb)      :: cloud_patch_threshold         !for the definition of cloud top
    
    !Grid data
    ! real(kind=kind_rb) :: total_cloud_fraction                    !total fraction of of grid that is covered by clouds
    ! real(kind=kind_rb),dimension(:,:) :: LWP_flattened     (imax, jmax)      !flattened collumns LWP content    
    ! real(kind=kind_rb),dimension(:,:,:) :: LWP_grid        (imax, jmax, krad1)    !full grid LWP contents, Is not actually a LWP
    ! real(kind=kind_rb),dimension(:,:) :: WVP_flattened     (imax, jmax)      !flattened collumns LWP content !Changed This (PIER_QV)
    !----------------------------------------------
    
    !Cloud ordering
    integer :: n1, n2                !Boundaries of start and end within GLQ
    integer :: n_classes                    !actual amount opf used classes, can be less then initial classes
    integer :: min_class              !amount of clouds in smallest cloud class
    integer :: counter,counter2            !counter that allows for cloud_top ordering
    real(kind=kind_rb) :: min_thresh              !Minimal size a cloud class has to have
    character(len=6) :: int_str_container                  !Is used to write ratio number into filenames

    
    real(kind=kind_rb),allocatable,dimension(:) :: quantiles_value           != n_classes - 1, could be integer but must be real for quicksort
    !This definition not necessary because it is defined later
    real(kind=kind_rb),allocatable,dimension(:) :: clear_WVP_ordered      !Ordered QV for cloudless collumns !Changed This (PIER_QV)
    real(kind=kind_rb),allocatable,dimension(:) :: cloudtop_height_ordered     !ordered cloudheights
    real(kind=kind_rb),allocatable,dimension(:,:) :: cloudtop_LWP_ordered    !Ordered LWP for cloudy collumns
    
    real(kind=kind_rb),allocatable,dimension(:) :: temp_quicksort_LWP_ordered
    integer,allocatable,dimension(:,:):: temp_quicksort_LWP_indexes
    
    real(kind=kind_rb),allocatable,dimension(:) :: temp_GLQ_points_cloudtop, temp_GLQ_weights_cloudtop
    
    real(kind=kind_rb),dimension(:,:,:) :: layerP_grid     (imax, jmax, krad1)    !pressure at grid core (full-level)
    real(kind=kind_rb),dimension(:,:,:) :: interfaceP_grid  (imax, jmax, krad2)    !pressure at grid interface (half-level)
    real(kind=kind_rb),dimension(:,:,:) :: layerMass_grid  (imax, jmax, krad1)    !mass within a gridpoint
    real(kind=kind_rb),dimension(:,:,:) :: qcl_grid      (imax, jmax, kradmax)  !actually just ql0
    
    
    real(kind=kind_rb),dimension(:,:,:) :: qv_grid      (imax, jmax, kradmax)                   !Changed This (PIER_QV)
    real(kind=kind_rb),dimension(:,:,:) :: WVP_grid      (imax, jmax, kradmax)                   !Changed This (PIER_QV)

    !This variable might nog be necessary
    integer,dimension(:) :: cloudtop_distribution (k1)           !Cloud height distribution amount of clouds in height index x
    real(kind=kind_rb),dimension(:,:) :: ztop_field  (imax, jmax)      !Grid containing collumn highest cloud, cloudtop height
    integer,dimension(:,:) :: cloud_class(imax, jmax)            !Contains the whole grid with integers showing to which class every collumns belongs
    
    real(KIND=kind_rb),dimension(:) :: LWP_distribution (kradmax)         !Column LWP distribution  
    real(kind=kind_rb),dimension(:,:) :: cloudFrac (imax, jmax)      ! cloud height class grid

    !__________________________________________________________
    !Define all field values
    
    qcl_grid(:, :, :) = 0.
    
    !Define qcl
    do k=1,kmax
      !im=i-1
      !jm=j-1
      qcl_grid(1:i1-1,1:j1-1,k) = ql0(2:i1,2:j1,k) !693
    end do

    do k=kmax+1,kradmax
      qcl_grid(:,:,k) = 0. !710
    end do
    
    !Define Pressure
    !Layer_P
    do k=1,kmax
      !im=i-1
      !jm=j-1
      layerP_grid(1:i1-1,1:j-1,k) = presf_input(k) !700
    end do
    
    do k=kmax+1,kradmax
      layerP_grid(1:imax,1:jmax,k)    = presf_input (k) !743
    end do
    layerP_grid(1:imax,1:jmax, krad1)   = 0.5*presh_input(krad1)!749
    
    !interface_P
    do k=1, krad1
      interfaceP_grid(1:imax,1:jmax,k) = presh_input(k) !764
    end do
    do i=1,imax
      interfaceP_grid(i,:, krad2)  = min( 1.e-4_kind_rb , 0.25*layerP_grid(1,:,krad1) ) !767
    end do

    !Define layerMass_grid
    layerMass_grid(1:imax,1:jmax,1:kradmax) = 100.*( interfaceP_grid(1:imax,1:jmax,1:kradmax) - interfaceP_grid(1:imax,1:jmax,2:kradmax+1) ) / grav  !of full level
    layerMass_grid(1:imax,1:jmax,krad1) = 100.*( interfaceP_grid(1:imax,1:jmax,krad1) - interfaceP_grid(1:imax,1:jmax,krad2) ) / grav
        
    !Define grid liquid water path
    LWP_grid(1:imax,1:jmax,1:kradmax) = qcl_grid(1:imax,1:jmax,1:kradmax)*layerMass_grid(1:imax,1:jmax,1:kradmax)*1e3
    LWP_grid(1:imax,1:jmax,krad1) = 0.

    !define qv and Water Vapor Path (WVP)
    do i=2,i1
      im=i-1
      do j=2,j1
        jm=j-1
        do k=1,kmax
          qv_grid(im,jm,k) = max(qt0(i,j,k) - ql0(i,j,k),1e-18) !avoid RRTMG reading negative initial values 
        enddo
        ksounding=npatch_start
        do k=kmax+1,kradmax
          qv_grid(im,jm,k) = qsnd(ksounding)
          ksounding=ksounding+1
        enddo
      enddo
    enddo

    WVP_grid(1:imax,1:jmax,1:kradmax) = qv_grid(1:imax,1:jmax,1:kradmax)*layerMass_grid(1:imax,1:jmax,1:kradmax)*1e3
    
    !Define the total LWP and WVP of every collumn
    !Define where the clouds are and at what height
    n_clear = 0
    n_clouds = 0
    cloudFrac(:,:)=0
    ztop_field(:,:) = 0
    cloudtop_distribution(:) = 0
    do i=1,imax
       do j=1,jmax
        ! Defines LWP_flattened and sets cloudFrac to 1 when LWP larger than cloud_threshold
        LWP_flattened(i,j) = SUM(LWP_grid(i,j,:))
        WVP_flattened(i,j) = SUM(WVP_grid(i,j,:))
        !Something could go wrong here with cloud threshold =/= cloud patch threshold,
        !there could be more n_clouds than SUM(cloudtop_distribution)
        if (LWP_flattened(i,j)>cloud_threshold) then
          cloudFrac(i,j) = 1
          n_clouds = n_clouds + 1
          do k=1,k1
            inverse_k = k1 + 1 - k
            !looks through the liquid in a gridpoint from top to bottom and assigns the first nonzero value to the ztop_field and then exits the vertical loop
            if (LWP_grid(i,j,inverse_k)>cloud_patch_threshold) then
              ztop_field(i,j) = zf(inverse_k) 
              cloudtop_distribution(inverse_k) = cloudtop_distribution(inverse_k)+1
              EXIT
            end if
          end do
        end if
       end do
      end do
    
    ! Write some lines if the amount of clouds in n_clouds does not coincide with the amount of cloudtops
    if (SUM(cloudtop_distribution).ne.n_clouds) then
      print *, "Warning: cloud patch threshold and cloud threshold have undeterminable results"
      print *, "SUM(cloudtop_distribution)"
      print *, SUM(cloudtop_distribution)
      print *, "n_clouds"
      print *, n_clouds
    end if
    
    !Define the amount of clear columns
    n_clear = (imax*jmax)-n_clouds
    
    total_cloud_fraction = float(n_clouds)/float(imax*jmax)
    
    !Perform the finding of GLQ points for the cloudless collumns
    temp_n_GLQ_clear = n_GLQ_clear
    if (n_clear > 0) then
      !Reduce the amount of GLQ points  if there are not enough clear collumns to house the GLQ points
      if (temp_n_GLQ_clear > n_clear) then
        temp_n_GLQ_clear = n_clear
      end if
      
      allocate (clear_WVP_ordered (n_clear))!Changed This (PIER_QV)
      allocate (original_clear_LWP_indexes (n_clear, 2))
      
      counter = 0
      !Place the original LWP values and actual coordinates into an array containing all the indexes.
      do j = 1, jmax
        do i = 1, imax
          if (LWP_flattened(i,j) <= cloud_threshold) then
            counter = counter + 1
            clear_WVP_ordered(counter) = WVP_flattened(i,j)
            !!Shift with a single index due to i=1 and j=1 being boundary values
            original_clear_LWP_indexes(counter, 1) = i + 1
            original_clear_LWP_indexes(counter, 2) = j + 1
          end if
        end do
      end do

      !Order on basis of WVP
      call quicksortindexes(clear_WVP_ordered, 1, n_clear, original_clear_LWP_indexes, n_clear)
    
      !Determine the indexes of the Gauss-Legendre points
      allocate (GLQ_points_clear   (temp_n_GLQ_clear))
      allocate (GLQ_weights_clear  (temp_n_GLQ_clear))
      allocate (GLQ_clear_LWP_indexes (temp_n_GLQ_clear, 2))
    
      !Determine GLQ points for the clear columns
      if (use_gauleg) then
        !Gauss-legendre Quadrature node choice
        call gauleg(float(1), float(n_clear), GLQ_points_clear, GLQ_weights_clear, temp_n_GLQ_clear)
      else
        !Evenly spaced in ordered WVP space
        if (use_evenly_spaced) then
        binwidth = float(n_clear)/float(temp_n_GLQ_clear)
          do nbin=1,temp_n_GLQ_clear
            GLQ_points_clear(nbin) = binwidth/2.0 + 0.5 + binwidth*(nbin-1)
          enddo
          GLQ_weights_clear(:) = 1.0 !Incorrect value, but not relevant
        else
          !GLQ points are determined on WVP bins
          if (use_bin) then
            valuewidth = (maxval(clear_WVP_ordered)-minval(clear_WVP_ordered))/float(temp_n_GLQ_clear)
            do nbin=1,temp_n_GLQ_clear
              GLQ_val = valuewidth/2.0 + minval(clear_WVP_ordered) + valuewidth*(nbin-1)
              temploc = minloc(abs(clear_WVP_ordered-GLQ_val))
              GLQ_points_clear(nbin) = temploc(1)
            enddo
            GLQ_weights_clear(:) = 1.0 !Incorrect value, but not relevant
          else
            !If nothing chosen, resort to GLQ selection
            call gauleg(float(1), float(n_clear), GLQ_points_clear, GLQ_weights_clear, temp_n_GLQ_clear)
          end if
        end if
      end if
    
      !Save coordinates of the points to an array containing all the clear column GLQ points
      do N_g = 1, temp_n_GLQ_clear
        x_index  = nint(GLQ_points_clear(N_g))
        temp_i   = int(original_clear_LWP_indexes(x_index,1))
        temp_j   = int(original_clear_LWP_indexes(x_index,2))
        GLQ_clear_LWP_indexes(N_g, 1) = temp_i
        GLQ_clear_LWP_indexes(N_g, 2) = temp_j
      end do
    
      deallocate (clear_WVP_ordered)
    else
      temp_n_GLQ_clear = 0
    end if

    !Cloudy Sky Gauss-Legendre
    !Select only the collumns with a nonzero cloudratio
    n_classes = 0
    n_GLQ_cloudtop = 0
    !Perform the finding of GLQ points for the clouded collumns
    if (n_clouds > 0) then
      allocate (cloudtop_height_ordered (n_clouds))
      
      !Place the clouds in order on basis of cloudtop height from low to high.
      counter = 0
      counter2 = 0
      do k=1,k1
        counter2 = cloudtop_distribution(k)
        do while (counter2>0)
          counter = counter + 1
          cloudtop_height_ordered(counter) = zf(k)
          counter2 = counter2 - 1
        end do
      end do
      
      !Initialise the classes
      !First tries to fill n_classes_initial with the same amount of collumns, 
      !if that fails it tries again with n_classes-1, this repeats until either all the classes are the same size or n_classes == 0
      !The classes are ordered on basis of cloudtop height
      n_classes = n_classes_initial
    10  if (n_classes > 1) then
          allocate (quantiles_value(n_classes-1))

        !Determines the edges for every class
        call quantiles (n_clouds, n_classes-1, .false., cloudtop_height_ordered, quantiles_value)

        allocate (n_in_class(n_classes))
        n_in_class(:) = 0
        cloud_class(:,:) = 0
      
        !determines the size of the classes and fills cloud_class with the integers of which the classes belong to
        do i = 1,imax
          do j = 1, jmax
            if(ztop_field(i,j) > 0) then
              if(ztop_field(i,j) > quantiles_value(n_classes-1)) then
                cloud_class(i,j) = n_classes
                n_in_class(n_classes) = n_in_class(n_classes) + 1
              else
                do n = 1, n_classes-1
                  if(ztop_field(i,j) <= quantiles_value(n)) then
                    cloud_class(i,j) = n
                    n_in_class(n) = n_in_class(n) + 1
                    exit
                  end if
                end do
              end if
            end if
          end do
        end do

        !If the "smallest" class is too small, retry making classes
        min_class  = minval(n_in_class(:))
        min_thresh = 0.01*float(imax*jmax) * total_cloud_fraction
        if (min_class < min_thresh) then    ! if too few in the least populated, reduce "n_classes" by 1 and redo...
          n_classes = n_classes - 1
          deallocate (quantiles_value)
          deallocate (n_in_class)
          goto 10
        end if
        
        !Checks if all the classes are the same size
        if (classes_same_size) then
          do i=2,n_classes
            if (n_in_class(1) /= n_in_class(i)) then
              n_classes = n_classes - 1
              deallocate (quantiles_value)
              deallocate (n_in_class)
              goto 10
            end if
          end do
        end if
        deallocate (quantiles_value)
      else
        allocate (n_in_class(n_classes))
        n_in_class(:) = 0
        n_in_class(1) = n_clouds
        cloud_class(:,:) = merge(1,0, cloudFrac>0)
      end if
      
      !Determine how many GLQ points have to be chosen for the cloudtop case
      n_RT = float(imax*jmax)/float(n_RT_Ratio)
      n_GLQ_cloudtop = nint(float(n_RT)/float(n_classes))
      if (n_GLQ_cloudtop <= 0) then
        n_GLQ_cloudtop = 1
      endif
      
      allocate(GLQ_in_class(n_classes))
      
      !Saves the amount of nodes(points) per class to {GLQ_in_class}
      if (dynamic_GLQ_per_class) then
        ! Amount of nodes per class dependent on amount of clouds per class
        if (min_GLQ_in_class>n_GLQ_cloudtop) then
          min_GLQ_in_class = n_GLQ_cloudtop
        end if
        do i=1,n_classes
          GLQ_in_class(i) = min_GLQ_in_class
          GLQ_in_class(i) = GLQ_in_class(i) + nint(float(n_in_class(i)*max((n_RT- min_GLQ_in_class*n_classes),0))/float(sum(n_in_class)))
        end do
        
        allocate(class_of_GLQ(sum(GLQ_in_class)))
        class_GLQ_counter = 1
        do i=1,sum(GLQ_in_class)
          if (i>sum(GLQ_in_class(1:class_GLQ_counter))) then
            class_GLQ_counter = class_GLQ_counter + 1
          end if
          class_of_GLQ(i) = class_GLQ_counter
        end do
      else
        ! Same amount of nodes in every class
        do i=1,n_classes
          GLQ_in_class(i) = n_GLQ_cloudtop
        end do
        allocate(class_of_GLQ(sum(GLQ_in_class)))
        class_GLQ_counter = 1
        do i=1,sum(GLQ_in_class)
          if (i>sum(GLQ_in_class(1:class_GLQ_counter))) THEN
            class_GLQ_counter = class_GLQ_counter + 1
          end if
          class_of_GLQ(i) = class_GLQ_counter
        end do
      end if
      
      allocate (GLQ_points_cloudtop (maxval(GLQ_in_class), n_classes))
      allocate (GLQ_weights_cloudtop(maxval(GLQ_in_class), n_classes))
      
      allocate (cloudtop_LWP_ordered(maxval(n_in_class), n_classes))
      allocate (original_cloudtop_LWP_indexes(maxval(n_in_class), 2, n_classes))
      
      allocate (GLQ_cloudtop_LWP_indexes(maxval(GLQ_in_class), 2, n_classes))
      
      do n = 1, n_classes
        !Determine Gauss-Legendre Quadrature points (nodes) for the clouded case      
        counter = 0
        do j = 1, jmax
          do i = 1, imax
            if (cloud_class(i,j) == n) then
              counter = counter + 1
              cloudtop_LWP_ordered(counter, n) = LWP_flattened(i, j)
              !!Shift with a single index due to i=1 and j=1 being boundary values
              original_cloudtop_LWP_indexes(counter, 1, n) = i + 1
              original_cloudtop_LWP_indexes(counter, 2, n) = j + 1
            end if
          end do
        end do
        !Sort the clouds on basis of LWP using quicksort, some other algorhitm could be used..
        !Save the indexes and values
        allocate(temp_quicksort_LWP_ordered(n_in_class(n)))
        allocate(temp_quicksort_LWP_indexes(n_in_class(n), 2))
        
        temp_quicksort_LWP_ordered = cloudtop_LWP_ordered(1:n_in_class(n),n)
        temp_quicksort_LWP_indexes = original_cloudtop_LWP_indexes(1:n_in_class(n),:,n)
        
        call quicksortindexes(temp_quicksort_LWP_ordered, 1, n_in_class(n), temp_quicksort_LWP_indexes, n_in_class(n))
        
        cloudtop_LWP_ordered(1:n_in_class(n),n) = temp_quicksort_LWP_ordered
        original_cloudtop_LWP_indexes(1:n_in_class(n),:,n) = temp_quicksort_LWP_indexes
        
        deallocate(temp_quicksort_LWP_ordered)
        deallocate(temp_quicksort_LWP_indexes)

        !Determine GLQ points for the clouded columns
        if (use_gauleg) then
          !Gauss-legendre Quadrature node choice
          allocate(temp_GLQ_points_cloudtop(GLQ_in_class(n)))
          allocate(temp_GLQ_weights_cloudtop(GLQ_in_class(n)))
          
          temp_GLQ_points_cloudtop = GLQ_points_cloudtop(1:GLQ_in_class(n), n)
          temp_GLQ_weights_cloudtop = GLQ_weights_cloudtop(1:GLQ_in_class(n), n)
          
          call gauleg(float(1), float(n_in_class(n)), temp_GLQ_points_cloudtop, temp_GLQ_weights_cloudtop, GLQ_in_class(n))
          
          GLQ_points_cloudtop(1:GLQ_in_class(n), n) = temp_GLQ_points_cloudtop
          GLQ_weights_cloudtop(1:GLQ_in_class(n), n) = temp_GLQ_weights_cloudtop
          
          deallocate(temp_GLQ_points_cloudtop)
          deallocate(temp_GLQ_weights_cloudtop)
        else
          !Evenly spaced in ordered WVP space
          if (use_evenly_spaced) then
            binwidth = float(n_in_class(n))/float(GLQ_in_class(n))
            do nbin=1,GLQ_in_class(n)
              GLQ_points_cloudtop(nbin, n) = binwidth/2.0 + 0.5 + binwidth*(nbin-1)
            enddo
            GLQ_weights_cloudtop(:, n) = 1.0 !Incorrect value, but not relevant
          else
            !GLQ points are determined on WVP bins
            if (use_bin) then
              valuewidth = (maxval(cloudtop_LWP_ordered(1:n_in_class(n),n))-minval(cloudtop_LWP_ordered(1:n_in_class(n),n)))/float(GLQ_in_class(n))
              do nbin=1,GLQ_in_class(n)
                GLQ_val = valuewidth/2.0 + minval(cloudtop_LWP_ordered(:,n)) + valuewidth*(nbin-1)
                temploc = minloc(abs(cloudtop_LWP_ordered(1:n_in_class(n),n)-GLQ_val))
                GLQ_points_cloudtop(nbin, n) = temploc(1)
              enddo
              GLQ_weights_cloudtop(:, n) = 1.0 !Incorrect value, but not relevant
            else
              !If nothing chosen, resort to GLQ selection
              allocate(temp_GLQ_points_cloudtop(GLQ_in_class(n)))
              allocate(temp_GLQ_weights_cloudtop(GLQ_in_class(n)))
              
              temp_GLQ_points_cloudtop = GLQ_points_cloudtop(1:GLQ_in_class(n), n)
              temp_GLQ_weights_cloudtop = GLQ_weights_cloudtop(1:GLQ_in_class(n), n)
              
              call gauleg(float(1), float(n_in_class(n)), temp_GLQ_points_cloudtop, temp_GLQ_weights_cloudtop, GLQ_in_class(n))
              
              GLQ_points_cloudtop(1:GLQ_in_class(n), n) = temp_GLQ_points_cloudtop
              GLQ_weights_cloudtop(1:GLQ_in_class(n), n) = temp_GLQ_weights_cloudtop
              
              deallocate(temp_GLQ_points_cloudtop)
              deallocate(temp_GLQ_weights_cloudtop)
            end if
          end if
        end if
        !Save coordinates of the points to an array containing all the clouded GLQ point indexes
        
        do N_g = 1, GLQ_in_class(n)  
          x_index = nint(GLQ_points_cloudtop(N_g, n))
          temp_i = int(original_cloudtop_LWP_indexes(x_index, 1, n))
          temp_j = int(original_cloudtop_LWP_indexes(x_index, 2, n))
          GLQ_cloudtop_LWP_indexes(N_g, 1, n) = temp_i
          GLQ_cloudtop_LWP_indexes(N_g, 2, n) = temp_j
        end do
      end do
      deallocate(cloudtop_height_ordered)
      deallocate(cloudtop_LWP_ordered)
    end if
    
    total_amount_GLQ_points = temp_n_GLQ_clear + sum(GLQ_in_class)    
    
    !!GLQ_indexes, indexes of the GLQ points
    allocate(GLQ_index_all(total_amount_GLQ_points, 2))

    !Places the clouded and clear GLQ points into a single array containing all the indexes of GLQ points
    if (temp_n_GLQ_clear>0) then
      do i =1, temp_n_GLQ_clear
        GLQ_index_all(i, 1) = GLQ_clear_LWP_indexes(i, 1)
        GLQ_index_all(i, 2) = GLQ_clear_LWP_indexes(i, 2)
      enddo
    end if
    GLQ_counter = temp_n_GLQ_clear

    do i=1,n_classes
      do j= 1,GLQ_in_class(i)
        GLQ_counter = GLQ_counter + 1
        GLQ_index_all(GLQ_counter, 1) = GLQ_cloudtop_LWP_indexes(j, 1, i)
        GLQ_index_all(GLQ_counter, 2) = GLQ_cloudtop_LWP_indexes(j, 2, i)
      enddo
    enddo

    if (n_clear >0) deallocate(GLQ_clear_LWP_indexes)
    if (n_clouds >0) deallocate(GLQ_cloudtop_LWP_indexes)

  end subroutine findGLQPoints
  
  !This subroutine places all the calculated values of GLQ points into the other points in the array
  ! subroutine reshuffleValues(n_GLQ_clear, GLQ_points_clear, GLQ_weights_clear, n_clear, &
    ! n_GLQ_cloudtop, GLQ_points_cloudtop, GLQ_weights_cloudtop, n_clouds, &
    ! n_classes,n_in_class, class_size, passed_GLQ_point, total_amount_GLQ_points, passed_slice_length, &
    ! original_clear_LWP_indexes, original_cloudtop_LWP_indexes)
  subroutine reshuffleValues(passed_GLQ_point, passed_slice_length)
  
    use modglobal, only: k1, boltz
    use modfields, only: thl0
      use modsurfdata, only: tskin
    use modraddata
    
    integer :: i, j
    integer :: fill_i, fill_j
    integer :: n, n1, n2
    integer :: class_number                              !Counter for the cloudtop classes
    integer :: passed_GLQ_point, temp_GLQ_point                  !GLQ point counter for the barker method
    integer :: cloudtop_GLQ_point                          !Reduced GLQ point counter for 
    integer :: passed_slice_length
    
    !!Contained in modraddata
    !----------------------------------------------
    ! integer :: n_GLQ_cloudtop, n_GLQ_clear                  !Amount of points for GLQ (cloudtop is per class)
    ! integer :: total_amount_GLQ_points                    !Total amount of GLQ points (n_GLQ_clear + n_GLQ_cloudtop*n_classes)
    ! integer,allocatable,dimension(:,:):: GLQ_index_all          !All GLQ indexes in a single array starting with cloudless and appending the first clouded class after being followed by second clouded etc.
    ! integer,allocatable,dimension(:,:):: original_clear_LWP_indexes       !original indexes of cloudless_LWP_ordered
    ! integer,allocatable,dimension(:,:,:):: original_cloudtop_LWP_indexes    !original indexes of the sorted LWP
    ! real(kind=kind_rb),allocatable,dimension(:) :: GLQ_points_clear, GLQ_weights_clear  !GLQ values cloudless
    ! real(kind=kind_rb),allocatable,dimension(:,:) :: GLQ_points_cloudtop, GLQ_weights_cloudtop  !GLQ values cloudtop, extra axis for the classes
    
    ! integer,allocatable,dimension(:) :: n_in_class                 !Array that contains the amount of clouds in a certain class
    ! integer :: class_size              !Amount of clouds in individual class
    ! integer :: n_clouds, n_clear          !number of collums with clouds and number of clear collumns
    ! integer :: n_classes_initial                !maximum number of cloudtop altitude classes
    ! real(kind=kind_rb)      :: cloud_threshold             !for the definition of a clouded collumn
    ! real(kind=kind_rb)      :: cloud_patch_threshold         !for the definition of cloud top
    
    !!This is for testpurposes, this makes it possible to check whether the program functions nicely
    ! integer,allocatable,dimension(:,:):: original_index_all
    !----------------------------------------------
  
    temp_GLQ_point = passed_GLQ_point
    !This loop passes through the slice and assigns values of a certain GLQ point to the pointss around that GLQ point as seen on basis of cloud height and collumn LWP
    do i=1,passed_slice_length
      
      if (temp_GLQ_point <= total_amount_GLQ_points) then
        if (temp_GLQ_point <= temp_n_GLQ_clear) then
          !Cloudless
          !Determine GLQ bin edges for replacing
          if (temp_n_GLQ_clear == 1) then
            n1=1
            n2=n_clear
          else
            if (temp_GLQ_point == 1) then
              n1 = 1
              n2 = nint((GLQ_points_clear(temp_GLQ_point) + GLQ_points_clear(temp_GLQ_point+1)) / 2)
            else
              if (temp_GLQ_point < temp_n_GLQ_clear) then
                n1 = nint((GLQ_points_clear(temp_GLQ_point-1) + GLQ_points_clear(temp_GLQ_point)) / 2)
                n1 = n1 + 1
                n2 = nint((GLQ_points_clear(temp_GLQ_point) + GLQ_points_clear(temp_GLQ_point+1)) / 2)
              else
                n1 = nint((GLQ_points_clear(temp_n_GLQ_clear-1) + GLQ_points_clear(temp_n_GLQ_clear)) / 2)
                n1 = n1 + 1
                n2 = n_clear
              end if
            end if
          end if
          !Places the values into the non GLQ nodes on the radiation arrays
          do n = n1, n2
            fill_i = int(original_clear_LWP_indexes(n, 1))
            fill_j = int(original_clear_LWP_indexes(n, 2))
            
            lwu(fill_i, fill_j,1:k1) =  lwUp_slice  (i,1:k1)
            lwd(fill_i, fill_j,1:k1) = -lwDown_slice(i,1:k1)
            !____________________
            if (.not. rad_longw) then !we get LW at surface identically to how it is done in sunray subroutine 
            !do i=2,i1
            lwd(fill_i, fill_j,1) =  -0.8 * boltz * thl0(fill_i, fill_j,1) ** 4.
            lwu(fill_i, fill_j,1) =  1.0 * boltz * tskin(fill_i, fill_j) ** 4.
            !end do
            end if
            !____________________
            swu(fill_i, fill_j,1:k1) =  swUp_slice  (i,1:k1)
            swd(fill_i, fill_j,1:k1) = -swDown_slice(i,1:k1)

            swdir(fill_i, fill_j,1:k1) = -swDownDir_slice(i,1:k1)
            swdif(fill_i, fill_j,1:k1) = -swDownDif_slice(i,1:k1)
            lwc  (fill_i, fill_j,1:k1) =  LWP_slice      (i,1:k1)

            lwuca(fill_i, fill_j,1:k1) =  lwUpCS_slice  (i,1:k1)
            lwdca(fill_i, fill_j,1:k1) = -lwDownCS_slice(i,1:k1)
            swuca(fill_i, fill_j,1:k1) =  swUpCS_slice  (i,1:k1)
            swdca(fill_i, fill_j,1:k1) = -swDownCS_slice(i,1:k1)

            SW_up_TOA (fill_i, fill_j) =  swUp_slice  (i,krad2)
            SW_dn_TOA (fill_i, fill_j) = -swDown_slice(i,krad2)
            
            LW_up_TOA (fill_i, fill_j) =  lwUp_slice  (i,krad2)
            LW_dn_TOA (fill_i, fill_j) = -lwDown_slice(i,krad2)

            SW_up_ca_TOA (fill_i, fill_j) =  swUpCS_slice  (i,krad2)
            SW_dn_ca_TOA (fill_i, fill_j) = -swDownCS_slice(i,krad2)
            
            LW_up_ca_TOA (fill_i, fill_j) =  lwUpCS_slice  (i,krad2)
            LW_dn_ca_TOA (fill_i, fill_j) = -lwDownCS_slice(i,krad2)
          end do
          temp_GLQ_point = temp_GLQ_point + 1
        else
          !cloudtop
          !Necessary to determine the clouded GLQ point with respect to the amount of clear GLQ points and to which class number the clouded GLQ point belongs 
          cloudtop_GLQ_point = temp_GLQ_point - temp_n_GLQ_clear
          
          class_number = class_of_GLQ(cloudtop_GLQ_point)
          do j=1,class_number-1
            cloudtop_GLQ_point = cloudtop_GLQ_point - GLQ_in_class(j)
          end do
              
          !Determine GLQ bin edges for replacing
          if (GLQ_in_class(class_number) == 1) then
            n1=1
            n2=n_in_class(class_number)
          else
            if (cloudtop_GLQ_point == 1) then
              n1 = 1
              n2 = nint((GLQ_points_cloudtop(cloudtop_GLQ_point, class_number) + GLQ_points_cloudtop(cloudtop_GLQ_point+1, class_number)) / 2)
            else
              if (cloudtop_GLQ_point < GLQ_in_class(class_number)) then
                n1 = nint((GLQ_points_cloudtop(cloudtop_GLQ_point-1, class_number) + GLQ_points_cloudtop(cloudtop_GLQ_point, class_number)) / 2)
                n1 = n1 + 1
                n2 = nint((GLQ_points_cloudtop(cloudtop_GLQ_point, class_number) + GLQ_points_cloudtop(cloudtop_GLQ_point+1, class_number)) / 2)
              else
                n1 = nint((GLQ_points_cloudtop(GLQ_in_class(class_number)-1, class_number) + GLQ_points_cloudtop(GLQ_in_class(class_number), class_number)) / 2)
                n1 = n1 + 1
                n2 = n_in_class(class_number)
              end if
            end if
          end if
          
          !Places the values into the non nodes places on the radiation arrays
          do n = n1, n2
            fill_i = int(original_cloudtop_LWP_indexes(n, 1, class_number))
            fill_j = int(original_cloudtop_LWP_indexes(n, 2, class_number))
            
            lwu(fill_i, fill_j,1:k1) =  lwUp_slice  (i,1:k1)
            lwd(fill_i, fill_j,1:k1) = -lwDown_slice(i,1:k1)

            !____________________
            if (.not. rad_longw) then !we get LW at surface identically to how it is done in sunray subroutine 
            !do i=2,i1
            lwd(fill_i, fill_j,1) =  -0.8 * boltz * thl0(fill_i, fill_j,1) ** 4.
            lwu(fill_i, fill_j,1) =  1.0 * boltz * tskin(fill_i, fill_j) ** 4.
            !end do
            end if
            !____________________

            swu(fill_i, fill_j,1:k1) =  swUp_slice  (i,1:k1)
            swd(fill_i, fill_j,1:k1) = -swDown_slice(i,1:k1)

            swdir(fill_i, fill_j,1:k1) = -swDownDir_slice(i,1:k1)
            swdif(fill_i, fill_j,1:k1) = -swDownDif_slice(i,1:k1)
            lwc  (fill_i, fill_j,1:k1) =  LWP_slice      (i,1:k1)

            lwuca(fill_i, fill_j,1:k1) =  lwUpCS_slice  (i,1:k1)
            lwdca(fill_i, fill_j,1:k1) = -lwDownCS_slice(i,1:k1)
            swuca(fill_i, fill_j,1:k1) =  swUpCS_slice  (i,1:k1)
            swdca(fill_i, fill_j,1:k1) = -swDownCS_slice(i,1:k1)

            SW_up_TOA (fill_i, fill_j) =  swUp_slice  (i,krad2)
            SW_dn_TOA (fill_i, fill_j) = -swDown_slice(i,krad2)
            LW_up_TOA (fill_i, fill_j) =  lwUp_slice  (i,krad2)
            LW_dn_TOA (fill_i, fill_j) = -lwDown_slice(i,krad2)

            SW_up_ca_TOA (fill_i, fill_j) =  swUpCS_slice  (i,krad2)
            SW_dn_ca_TOA (fill_i, fill_j) = -swDownCS_slice(i,krad2)
            LW_up_ca_TOA (fill_i, fill_j) =  lwUpCS_slice  (i,krad2)
            LW_dn_ca_TOA (fill_i, fill_j) = -lwDownCS_slice(i,krad2)
          end do
          temp_GLQ_point = temp_GLQ_point + 1
        end if
      end if
    end do
  end subroutine reshuffleValues
  
end module modradrrtmgyuri