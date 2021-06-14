module modradrrtmgyuri
use modraddata

implicit none

contains
!__________________________________________________________
!__________________________________________________________
!__________________________________________________________
!__________________________________________________________
!Important Routine
	subroutine findGLQPoints(n_GLQ_clear, GLQ_points_clear, GLQ_weights_clear, GLQ_clear_LWP_indexes, n_clear, &
		n_GLQ_cloudtop, GLQ_points_cloudtop, GLQ_weights_cloudtop, GLQ_cloudtop_LWP_indexes, n_clouds, &
		n_classes, n_class, class_size, total_amount_GLQ_points, GLQ_index_all, &
		original_clear_LWP_indexes, original_cloudtop_LWP_indexes)
	
		use modglobal, only: imax, jmax, kmax, i1, j1, k1, kind_rb, grav, zf, zh
		use modfields, only: ql0
		use modradrrtmgyuriroutines
		
		implicit none
		
		integer :: i,j,k,inverse_k, n
		
		!acceleration
		integer :: n_RT_Ratio 						!Ratio of Radiative transfer acceleration
		integer :: n_RT != nx * ny / n_RT_Ratio		!Actual number of Radtiative transfer calculations
		
		!Gauss-Legendre Quadrature
		integer :: N_g															!Loops through GLQ_points
		integer :: x_index														!Index of GLQ x_i
		integer :: temp_i, temp_j												!Temporary placeholders for i and j in sorting and ordering
		integer :: n_GLQ_cloudtop, n_GLQ_clear									!Amount of points for GLQ (cloudtop is per class)
		integer :: GLQ_counter													!counter necessary for full list of GLQ indexes
		integer :: total_amount_GLQ_points										!Total amount of GLQ points (n_GLQ_clear + n_GLQ_cloudtop*n_classes)
		real(kind=kind_rb),allocatable,dimension(:) :: GLQ_points_clear, GLQ_weights_clear	!GLQ values cloudless
		real(kind=kind_rb),allocatable,dimension(:,:) :: GLQ_points_cloudtop, GLQ_weights_cloudtop	!GLQ values cloudtop, extra axis for the classes
		integer,allocatable,dimension(:,:):: GLQ_clear_LWP_indexes  		!The indexes necessary for the GLQ of the cloudless LWP
		integer,allocatable,dimension(:,:):: GLQ_index_all					!All GLQ indexes in a single array starting with cloudless and appending the first clouded class after being followed by second clouded etc.
		integer,allocatable,dimension(:,:,:):: GLQ_cloudtop_LWP_indexes   		!The indexes necessary for the GLQ of the cloudtop LWP, extra axis for the classes
		integer,allocatable,dimension(:,:):: original_clear_LWP_indexes 			!original indexes of cloudless_LWP_ordered
		integer,allocatable,dimension(:,:,:):: original_cloudtop_LWP_indexes		!original indexes of the sorted LWP
		!integer,allocatable,dimension(:,:):: GLQ_cloudtop_height_indexes
		!integer,allocatable,dimension(:,:):: original_cloudtop_height_indexes		!original indexes of cloudtop_LWP_ordered
		
		!Cloud ordering
		integer :: n1, n2								!Boundaries of start and end within GLQ
		integer :: class_size							!Amount of clouds in individual class
		integer :: n_clouds, n_clear					!number of collums with clouds and number of clear collumns
		integer :: n_classes		            		!actual amount opf used classes, can be less then initial classes
		integer :: n_classes_initial            		!maximum number of cloudtop altitude classes
		integer :: min_class							!amount of clouds in smallest cloud class
		integer :: counter,counter2						!counter that allows for cloud_top ordering
		real(kind=kind_rb)    	:: cloud_threshold 						!for the definition of a clouded collumn
		real(kind=kind_rb)    	:: cloud_patch_threshold 				!for the definition of cloud top
		real(kind=kind_rb)		:: min_thresh							!Minimal size a cloud class has to have
		
		!Cloud ordering
		!integer,allocatable,dimension(:) :: clear_LWP_distribution					!Possibly very unnecessary
		!integer,allocatable,dimension(:) :: LWP_distribution						!Possibly very unneccesary
		integer,allocatable,dimension(:) :: n_class 								!Array that contains the amount of clouds in a certain class
		
		real(kind=kind_rb),allocatable,dimension(:) :: quantiles_value 					!= n_classes - 1, could be integer but must be real for quicksort
		!This definition not necessary because it is defined later
		!real(kind=kind_rb),allocatable,dimension(:) :: cloudtop_distribution 		!amount of cloudtops in every row, could be integer but must be real for quicksort
		real(kind=kind_rb),allocatable,dimension(:) :: clear_LWP_ordered			!Ordered LWP for cloudless collumns
		real(kind=kind_rb),allocatable,dimension(:) :: cloudtop_height_ordered 		!ordered cloudheights
		real(kind=kind_rb),allocatable,dimension(:,:) :: cloudtop_LWP_ordered		!Ordered LWP for cloudy collumns
		!real(kind=kind_rb),allocatable,dimension(:,:) :: subset_holder !!???
		
		!Grid data
		real(kind=kind_rb) :: total_cloud_fraction										!total fraction of of grid that is covered by clouds
		real(kind=kind_rb),dimension(:,:) :: LWP_flattened 		(imax, jmax)			!flattened collumns LWP content
		real(kind=kind_rb),dimension(:,:,:) :: LWP_grid    		(imax, jmax, krad1)		!full grid LWP contents, Is not actually a LWP
		real(kind=kind_rb),dimension(:,:,:) :: layerP_grid 		(imax, jmax, krad1)		!pressure at grid core (full-level)
		real(kind=kind_rb),dimension(:,:,:) :: interfaceP_grid	(imax, jmax, krad2)		!pressure at grid interface (half-level)
		real(kind=kind_rb),dimension(:,:,:) :: layerMass_grid	(imax, jmax, krad1)		!mass within a gridpoint
		real(kind=kind_rb),dimension(:,:,:) :: qcl_grid			(imax, jmax, kradmax)	!actually just ql0

		!This variable might nog be necessary
		real(KIND=kind_rb),dimension(:) :: cloudtop_distribution (k1)           !Cloud height distribution amount of clouds in height index x
		real(kind=kind_rb),dimension(:,:) :: ztop_field	(imax, jmax)			!Grid containing collumn highest cloud, cloudtop height
		integer,dimension(:,:) :: cloud_class(imax, jmax)						!Contains the whole grid with integers showing to which class every collumns belongs

		real(KIND=kind_rb),dimension(:) :: LWP_distribution (kradmax)         !Cloud height distribution	
		real(kind=kind_rb),dimension(:,:) :: cloudFrac (imax, jmax)			! cloud height class grid

		!__________________________________________________________
		!Define all field values
		!print *, "Starting Define all field values"
		!Should initialize all the variables
		!Make checklist
		
		!Things that are/should already be known should be defined here:
		!For example the allocation of space for arrays
		!Also, 
		n_RT_Ratio = 1
		n_classes_initial = 20
		cloud_threshold = 0.0
		cloud_patch_threshold = 0.0
		
		
		qcl_grid(:, :, :) = 0.
		!Maybe shift all instances of gl0
		
		!Deze sequence kan efficienter maar voor nu zoveel mogelijk gescheiden om implosies te voorkomen
		
		!Moet j ook nog verschoven worden of niet? Want hier laat ik j lopen van 1 tot jmax, maar is dat de bedoeling?
		
		!Checken of ik nog wat dingen moet opschuiven
		!Define qcl
		!do i=2,i1
		!do j=2,j1
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
		!do i=2,i1
		!do j=2,j1
		do k=1,kmax
			!im=i-1
			!jm=j-1
			layerP_grid(1:i1-1,1:j-1,k) = presf_input(k) !700
		end do
		
		!do i=1,imax
		!do j = 1,jmax
		do k=kmax+1,kradmax
			layerP_grid(1:imax,1:jmax,k)    = presf_input (k) !743
		end do
		layerP_grid(1:imax,1:jmax, krad1)   = 0.5*presh_input(krad1)!749

		
		!interface_P
		!do i=1,imax
		!do j = 1,jmax
		do k=1, krad1
			interfaceP_grid(1:imax,1:jmax,k) = presh_input(k) !764
		end do
		!end do
		!!!!Hier dus die : en die wil ik vervangen door iets explicieters, 1:j1 bijv.
		do i=1,imax
			interfaceP_grid(i,:, krad2)  = min( 1.e-4_kind_rb , 0.25*layerP_grid(1,:,krad1) ) !767
		end do
		!end do
			
		!do i=1,imax
	    !   interfaceP_grid(i,:, krad2)  = min( 1.e-4_kind_rb , 0.25*layerP_grid(1,:,krad1) )
		!end do
	    
	    !do k=kmax+1, kradmax
		!   layerP_grid(:,:,k) = presf_input(k)
	    !end do
	    !layerP_grid (:,:, krad1) = 0.5*presh_input(krad1)
	
		!This one could be partially moved to the other k do loops
		!Define Layermass and LWP
		
		!!!!Checken of er niet over i en j gesommeerd moet worden in plaats van aannemen dat het geheel wel klopt

		!777
		layerMass_grid(1:imax,1:jmax,1:kradmax) = 100.*( interfaceP_grid(1:imax,1:jmax,1:kradmax) - interfaceP_grid(1:imax,1:jmax,2:kradmax+1) ) / grav  !of full level
		!778
		LWP_grid(1:imax,1:jmax,1:kradmax) = qcl_grid(1:imax,1:jmax,1:kradmax)*layerMass_grid(1:imax,1:jmax,1:kradmax)*1e3

		!782 (possibly unneccesary)
	    layerMass_grid(1:imax,1:jmax,krad1) = 100.*( interfaceP_grid(1:imax,1:jmax,krad1) - interfaceP_grid(1:imax,1:jmax,krad2) ) / grav
        !784
		LWP_grid(1:imax,1:jmax,krad1) = 0.
		
		call writetofiledefinedsize("qcl_grid", qcl_grid, 3, imax, jmax, kradmax)
		call writetofiledefinedsize("layerP_grid", layerP_grid, 3, imax, jmax, krad1)
		call writetofiledefinedsize("interfaceP_grid", interfaceP_grid, 3, imax, jmax, krad2)
		call writetofiledefinedsize("layerMass_grid", layerMass_grid, 3, imax, jmax, krad1)
		call writetofiledefinedsize("LWP_grid", LWP_grid, 3, imax, jmax, krad1)
		
		!print *, "finished Define all field values"
		!__________________________________________________________
		!determine specific necessary values
		!print *, "starting cloud and LWP data"
		! cloud_threshold = 0.0
		! cloud_patch_threshold = 0.0
		
		
	    !LWP_grid is not cumulative I think, Also I think it has different values for different z
	    !I need to flatten it.
		!Define Ztop and cloudfrac
		n_clear = 0
		n_clouds = 0
		cloudFrac(:,:)=0
	    ztop_field(:,:) = 0
		cloudtop_distribution(:) = 0
		do i=1,imax
		   do j=1,jmax
				! Defines LWP_flattened and sets cloudFrac to 1 when LWP larger than cloud_threshold
				LWP_flattened(i,j) = SUM(LWP_grid(i,j,:))
				!Something could go wrong here with cloud threshold =/= cloud patch threshold,
				!there could be more n_clouds than SUM(cloudtop_distribution)
				if (LWP_flattened(i,j)>cloud_threshold) then
					cloudFrac(i,j) = 1
					n_clouds = n_clouds + 1
					do k=1,k1
						inverse_k = k1 + 1 - k
						!looks through the liquid in a gridpoint from top to bottom and assigns the first nonzero value to the ztop_field and then exits the vertical loop
						!if (ql0(i,j,inverse_k)>cloud_patch_threshold) then
						if (LWP_grid(i,j,inverse_k)>cloud_patch_threshold) then
							ztop_field(i,j) = zf(inverse_k) 
							cloudtop_distribution(inverse_k) = cloudtop_distribution(inverse_k)+1.0
							EXIT
						end if
					end do
				end if
		   end do
	    end do
		
		if (SUM(cloudtop_distribution).ne.n_clouds) then
			! print *, "Warning: cloud patch threshold and cloud threshold have undeterminable results"
			! print *, "SUM(cloudtop_distribution)"
			! print *, SUM(cloudtop_distribution)
			! print *, "n_clouds"
			! print *, n_clouds
		end if
		
		n_clear = (imax*jmax)-n_clouds
		
		total_cloud_fraction = float(n_clouds)/float(imax*jmax)
		
		
		call writetofiledefinedsize("ztop_field", ztop_field, 2, imax, jmax, 1)
		call writetofiledefinedsize("cloudtop_distribution", cloudtop_distribution, 1, k1, 1, 1)
		
		!!Cloudtop distribution might be unccecessary
		!!!~The sorting is bad, also you should not pass integers to quicksort!
		! print *, "cloudtop_distribution"
		! print *, cloudtop_distribution
		!call quicksort(cloudtop_distribution, 1, k1, k1)
		! print *, "sorted cloudtop_distribution"
		! print *, cloudtop_distribution
		
		! print *, "n_clouds"
		! print *, n_clouds
		! print *, "n_clear"
		! print *, n_clear
		! print *, "total_cloud_fraction"
		! print *, total_cloud_fraction
		
		!Determined:
		!   Cloud_fraction
		!   ztop_field
		!   Cloudtop_distibution
		!   n_clouds
		!   n_clear
		!   LWP_flattened
		
		! print *, "finished cloud and LWP data"
		! print *, "starting cloudless collumns"
		
		n_GLQ_clear = 0
		

		if (n_clear > 0) then
			n_GLQ_clear = 30
			if (n_GLQ_clear < n_clear) then
				n_GLQ_clear = n_clear
			end if
			! print *, "n_clear > 0"
			allocate (clear_LWP_ordered (n_clear))
			allocate (original_clear_LWP_indexes (n_clear, 2))
			
			counter = 0
			!!This could be moved up
			do i = 1, imax
				do j = 1, jmax
					if (LWP_flattened(i,j) <= cloud_threshold) then
						counter = counter + 1
						clear_LWP_ordered(counter) = LWP_flattened(i,j)
						!!Shift with a single index due to i=1 and j=1 being boundary values
						original_clear_LWP_indexes(counter, 1) = i + 1
						original_clear_LWP_indexes(counter, 2) = j + 1
						! original_clear_LWP_indexes(counter, 1) = i
						! original_clear_LWP_indexes(counter, 2) = j
					end if
				end do
			end do

			!Order on basis of LWP
			call quicksortindexes(clear_LWP_ordered, 1, n_clear, original_clear_LWP_indexes, n_clear)
		
			!Determine the indexes of the Gauss-Legendre points
			allocate (GLQ_points_clear 	(n_GLQ_clear))
			allocate (GLQ_weights_clear	(n_GLQ_clear))
			allocate (GLQ_clear_LWP_indexes (n_GLQ_clear, 2))
		
			call gauleg(float(1), float(n_clear), GLQ_points_clear, GLQ_weights_clear, n_GLQ_clear)
		
			do N_g = 1, n_GLQ_clear
			  x_index  = nint(GLQ_points_clear(N_g))
			  temp_i   = int(original_clear_LWP_indexes(x_index,1))
			  temp_j   = int(original_clear_LWP_indexes(x_index,2))
			  GLQ_clear_LWP_indexes(N_g, 1) = temp_i
			  GLQ_clear_LWP_indexes(N_g, 2) = temp_j
			  
			  !!NEED THE FLUX HERE
			  !F_GLQ = F_full(temp_i, temp_j)
			end do
		
			!deallocate (GLQ_points_clear)
			!deallocate (GLQ_weights_clear)
		end if
		! print *, "finished cloudless collumns"
		!____________________!!!!!!!!!!!!!!!!!_____________________
		!Cloudy Sky Gauss-Legendre
		!____________________!!!!!!!!!!!!!!!!!_____________________
		!Select only the collumns with a nonzero cloudratio
		! print *, "starting clouded collumns"
		if (n_clouds > 0) then
			! print *, "n_clouds > 0"
			allocate (cloudtop_height_ordered (n_clouds))

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
			
			call writetofiledefinedsize("LookIfActuallyIncreasing", cloudtop_height_ordered, 1, n_clouds, 1, 1)
			
			!Determined:
			!   cloudtop_height_ordered

			! print *, "cloudtop_height_ordered"
			! print *, cloudtop_height_ordered
			
			write(*,*)

			!Initialise the classes
			
			!allocate(cloud_class(imax, jmax))

			n_classes = n_classes_initial
		10	if (n_classes > 1) then
		    	allocate (quantiles_value(n_classes-1))

				! write(*,*) 'n_classes = ', n_classes

				call quantiles (n_clouds, n_classes-1, .false., cloudtop_height_ordered, quantiles_value)

				! print *, "quantiles_value"
				! print *, quantiles_value

				allocate (n_class(n_classes))
				n_class(:) = 0
				cloud_class(:,:) = 0
			
				do i = 1,imax
					do j = 1, jmax
						if(ztop_field(i,j) > 0) then
							if(ztop_field(i,j) > quantiles_value(n_classes-1)) then
								cloud_class(i,j) = n_classes
								n_class(n_classes) = n_class(n_classes) + 1
							else
								do n = 1, n_classes-1
									if(ztop_field(i,j) <= quantiles_value(n)) then
										cloud_class(i,j) = n
										n_class(n) = n_class(n) + 1
									end if
								end do
							end if
						end if
					end do
				end do
				
				!counter = 1
				!do n = 1,n_classes-1
				!	do while(cloudtop_height_ordered(counter) <= quantiles_value(n))
				!		n_class(n) = n_class(n) + 1
				!		counter = counter + 1
				!	end do
				!end do
				!counter = counter - 1
				!n_class(n_classes) = n_clouds - counter

				min_class  = minval(n_class(:))
				!Skip this part for now, but might need to look at this later
				min_thresh = 0.01*float(imax*jmax) * total_cloud_fraction
				if (min_class < min_thresh) then    ! if too few in the least populated, reduce "n_classes" by 1 and redo...
					n_classes = n_classes - 1
					deallocate (quantiles_value)
					deallocate (n_class)
					goto 10
				end if
				
				!print *, ztop_field(:,:)


				! write(*,*)
				! write(*,*) 'altitude of classes:'
				! do i = 1, n_classes
					! if (i < n_classes) then
						! write(*,*) "i, quantiles_value(i)"
						! write(*,*) i, quantiles_value(i)
						! write(*,*) "i, n_class(i)"
						! write(*,*) i, n_class(i)
					! else
						! write(*,*) "i, maxval(cloudtop_height_ordered)"
						! write(*,*) i, maxval(cloudtop_height_ordered)
						! write(*,*) "i, n_class(i)"
						! write(*,*) i, maxval(cloudtop_height_ordered)
					! end if
				! end do
				! write(*,*)
				
				class_size = n_class(1)
				do i=2,n_classes
					if (class_size /= n_class(i)) then
						print *, "WARNING: Something went wrong with cloud allocation in modradrrtmg, reducing class size"
						!!Maybe not necessary to stop loop?
						n_classes = n_classes - 1
						deallocate (quantiles_value)
						deallocate (n_class)
						goto 10
					end if
				end do
				!deallocate (quantiles_value)
			else
				allocate (n_class(n_classes))
				n_class(:) = 0
				n_class(1) = n_clouds
				class_size = n_class(1)
				cloud_class(:,:) = n_classes
			end if
			! print *, "cloud and classes finished"
			

			n_RT = (imax*jmax)/(n_RT_Ratio)
			n_GLQ_cloudtop = nint(float(n_RT)/float(n_classes))
			
			! print *, "GLQ_points_cloudtop"
			allocate (GLQ_points_cloudtop (n_GLQ_cloudtop, n_classes))
			! print *, "GLQ_weight_cloudtop"
			allocate (GLQ_weights_cloudtop(n_GLQ_cloudtop, n_classes))
			
			! print *, "cloudtop_LWP_ordered"
			allocate (cloudtop_LWP_ordered(class_size, n_classes))
			! print *, "original_cloudtop_LWP_indexes"
			allocate (original_cloudtop_LWP_indexes(class_size, 2, n_classes))
			
			! print *, "GLQ_cloudtop_LWP_indexes"
			allocate (GLQ_cloudtop_LWP_indexes(n_GLQ_cloudtop, 2, n_classes))
			
			! print *, "start going through classes"
			do n = 1, n_classes
				! print *, "n"
				! print *, n
				
				! print *, "gauleg"
				call gauleg(float(1), float(n_class(n)), GLQ_points_cloudtop(:, n), GLQ_weights_cloudtop(:, n), n_GLQ_cloudtop)
				
				
				
				! print *, "placing in original indexes"
				counter = 0
				do i = 1, imax
					do j = 1, jmax
						if (cloud_class(i,j) == n) then
							counter = counter + 1
							cloudtop_LWP_ordered(counter, n) = LWP_flattened(i, j)
							!!Shift with a single index due to i=1 and j=1 being boundary values
							! print *, counter
							! print *, i + 1
							! print *, j + 1
							original_cloudtop_LWP_indexes(counter, 1, n) = i + 1
							original_cloudtop_LWP_indexes(counter, 2, n) = j + 1
							! original_cloudtop_LWP_indexes(counter, 1, n) = i
							! original_cloudtop_LWP_indexes(counter, 2, n) = j
						end if
					end do
				end do
				
				call writetofiledefinedsize("cloudtop_LWP_ordered", cloudtop_LWP_ordered, 1, class_size, 1, 1)
				call writetofiledefinedsize("original_cloudtop_LWP_indexes", original_cloudtop_LWP_indexes, 1, class_size, 2, 1)
				
				! print *, "quicksortindexes"
				!Removed this for tests
				!!call quicksortindexes(cloudtop_LWP_ordered(:,n), 1, class_size, original_cloudtop_LWP_indexes(:,:,n), class_size)
			
				! print *, "save GLQ points"
				n2 = 0
				do N_g = 1, n_GLQ_cloudtop
				
					!!TODO: NOTHING DONE WITH n here
					! n1 = n2 + 1
					! if (N_g < n_GLQ_cloudtop) then
						! n2 = (GLQ_points_cloudtop(N_g, n) + GLQ_points_cloudtop(N_g + 1, n)) / 2
					! else
						! n2 = n_class(n)
					! end if
					
					

					!!print *, "bars are set, now placing  in GLQ indexes"
					!Look at if this works, weird index results
					!!!x_index = nint(GLQ_points_cloudtop(N_g, n))
					x_index = N_g
					temp_i = int(original_cloudtop_LWP_indexes(x_index, 1, n))
					temp_j = int(original_cloudtop_LWP_indexes(x_index, 2, n))
					GLQ_cloudtop_LWP_indexes(N_g, 1, n) = temp_i
					GLQ_cloudtop_LWP_indexes(N_g, 2, n) = temp_j
					
				end do
				
			end do

		end if
		! print *, "finished clouded collumns"
		! print *, "original_cloudtop_LWP_indexes(:,:,:)"
		! print *, original_cloudtop_LWP_indexes(:,:,:)
		!!!It might be unneccesary to make a total thing... ///  https://michaelgoerz.net/notes/advanced-array-passing-in-fortran.html
		! print *, "starting GLQ to long total array"
		total_amount_GLQ_points = n_GLQ_clear + n_GLQ_cloudtop*n_classes
		
		! print *, "allocating this amount of points", total_amount_GLQ_points
		allocate(GLQ_index_all(total_amount_GLQ_points, 2))
		
		! print *, "GLQ clear"
		if (n_GLQ_clear>0) then
			do i =1, n_GLQ_clear
				GLQ_index_all(i, 1) = GLQ_clear_LWP_indexes(i, 1)
				GLQ_index_all(i, 2) = GLQ_clear_LWP_indexes(i, 2)
			enddo
		end if
		GLQ_counter = n_GLQ_clear
		! print *, "GLQ clouded"
		do i=1,n_classes
			do j= 1,n_GLQ_cloudtop
				GLQ_counter = GLQ_counter + 1
				GLQ_index_all(GLQ_counter, 1) = GLQ_cloudtop_LWP_indexes(j, 1, i)
				GLQ_index_all(GLQ_counter, 2) = GLQ_cloudtop_LWP_indexes(j, 2, i)
			enddo
		enddo

		!!!print *, GLQ_points_cloudtop(:, 1)
		!print *, GLQ_index_all(:, :)
		
		! print *, "finished GLQ to long total array"
		!!!!!RADIATION!!!!!
		
		!__________________________________________________________
		!____________________!!!!!!!!!!!!!!!!!_____________________
		!FILL ARRAYS WITH SUBSETS
		!____________________!!!!!!!!!!!!!!!!!_____________________
		!__________________________________________________________
		
		
		!__________________________________________________________
		!Make and write to files
		!call writetofiles("ql0content", ql0, 3)
		
	end subroutine findGLQPoints
	
	!Only the values in the radiation or also the field values in modraddata?
	subroutine reshuffleValues(n_GLQ_clear, GLQ_points_clear, GLQ_weights_clear, GLQ_clear_LWP_indexes, n_clear, &
		n_GLQ_cloudtop, GLQ_points_cloudtop, GLQ_weights_cloudtop, GLQ_cloudtop_LWP_indexes, n_clouds, &
		n_classes,n_class, class_size, passed_GLQ_point, total_amount_GLQ_points, passed_slice_length, &
		original_clear_LWP_indexes, original_cloudtop_LWP_indexes)
	
	use modglobal, only: k1, boltz
	use modfields, only: thl0
    use modsurfdata, only: tskin
	use modraddata
	
	integer :: i
	!integer :: GLQ_i, GLQ_j
	integer :: fill_i, fill_j
	integer :: n, n1, n2
	integer :: class_number															!Counter for the cloudtop classes
	integer :: passed_GLQ_point, temp_GLQ_point									!GLQ point counter for the barker method
	integer :: cloudtop_GLQ_point													!Reduced GLQ point counter for 
	integer :: n_GLQ_cloudtop, n_GLQ_clear											!Amount of points for GLQ
	integer :: total_amount_GLQ_points										!Total amount of GLQ points (n_GLQ_clear + n_GLQ_cloudtop*n_classes)
	integer :: passed_slice_length
	real(kind=kind_rb),allocatable,dimension(:) :: GLQ_points_clear, GLQ_weights_clear	!GLQ values cloudless
	real(kind=kind_rb),allocatable,dimension(:,:) :: GLQ_points_cloudtop, GLQ_weights_cloudtop	!GLQ values cloudtop, extra axis for the classes
	integer,allocatable,dimension(:,:):: GLQ_clear_LWP_indexes  		!The indexes necessary for the GLQ of the cloudless LWP
	!integer,allocatable,dimension(:,:):: GLQ_index_all					!All GLQ indexes in a single array starting with cloudless and appending the first clouded class after being followed by second clouded etc.
	integer,allocatable,dimension(:,:,:):: GLQ_cloudtop_LWP_indexes   		!The indexes necessary for the GLQ of the cloudtop LWP, extra axis for the classes
	integer,allocatable,dimension(:,:):: original_clear_LWP_indexes 			!original indexes of cloudless_LWP_ordered
	integer,allocatable,dimension(:,:,:):: original_cloudtop_LWP_indexes		!original indexes of cloudtop_LWP_ordered


	integer,allocatable,dimension(:) :: n_class 								!Array that contains the amount of clouds in a certain class
	integer :: class_size							!Amount of clouds in individual class
	integer :: n_clouds, n_clear					!number of collums with clouds and number of clear collumns
	integer :: n_classes		            		!actual amount of used classes, can be less then initial classes
	integer :: n_classes_initial            		!maximum number of cloudtop altitude classes
	real    :: cloud_threshold 						!for the definition of a clouded collumn
	real    :: cloud_patch_threshold 				!for the definition of cloud top
		
		temp_GLQ_point = passed_GLQ_point
		do i=1,passed_slice_length
			
			!Is GLQ_index_all even necessary?
			!GLQ_i   = int(GLQ_index_all(temp_GLQ_point, 1))
			!GLQ_j   = int(GLQ_index_all(temp_GLQ_point, 2))
			
			! N_g = MODULO(temp_GLQ_point, passed_slice_length) + 1 !This is just i?
			if (temp_GLQ_point <= total_amount_GLQ_points) then
				! print *, "temp_GLQ_point <= total_amount_GLQ_points"
				if (temp_GLQ_point <= n_GLQ_clear) then
					! print *, "temp_GLQ_point <= n_GLQ_clear"
				!!!!!!!!!!!!!!!!!!!!!
					!Cloudless
					!n1 and n2 could be saved..., so you dont have to redetermine the n1 and n2
					if (temp_GLQ_point == 1) then
						n1 = 1
						n2 = (GLQ_points_clear(temp_GLQ_point) + GLQ_points_clear(temp_GLQ_point+1)) / 2
					else
							!!TODO: n_GLQ_clear is fine for now but should change when clear is put into classes
						if (temp_GLQ_point < n_GLQ_clear) then
							n1 = (GLQ_points_clear(temp_GLQ_point-1) + GLQ_points_clear(temp_GLQ_point)) / 2
							n2 = (GLQ_points_clear(temp_GLQ_point) + GLQ_points_clear(temp_GLQ_point+1)) / 2
						else
							! print *, "this is not happening right?"
							n1 = (GLQ_points_clear(n_clear-1) + GLQ_points_clear(n_clear)) / 2
							n2 = n_clear
						end if
					end if

					!F_GLQ = flux(temp_i, temp_j)
					do n = n1, n2
						!Is looking at original Clear LWP indexes correct? Have they been ordered?
						fill_i = int(original_clear_LWP_indexes(n, 1))
						fill_j = int(original_clear_LWP_indexes(n, 2))
						
						lwu(fill_i, fill_j,1:k1) =  lwUp_slice  (i,1:k1)
						lwd(fill_i, fill_j,1:k1) = -lwDown_slice(i,1:k1)
						!!How do this?
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
					! print *, "temp_GLQ_point > n_GLQ_clear"
					!cloudtop

					cloudtop_GLQ_point = temp_GLQ_point - n_GLQ_clear
					!!!
					class_number = cloudtop_GLQ_point/class_size
					if (MODULO(cloudtop_GLQ_point, class_size) > 0) then
						class_number = class_number + 1
					end if
					cloudtop_GLQ_point = cloudtop_GLQ_point - (class_number-1)*class_size
					!!cloudtop_GLQ_point = cloudtop_GLQ_point + 1
					
					! print *, "class_size"
					! print *, class_size
					! print *, "class_number"
					! print *, class_number
					! print *, "cloudtop_GLQ_point"
					! print *, cloudtop_GLQ_point
					! if (cloudtop_GLQ_point == 1) then
						! n1 = 1
						! n2 = (GLQ_points_cloudtop(cloudtop_GLQ_point, class_number) + GLQ_points_cloudtop(cloudtop_GLQ_point+1, class_number)) / 2
					! else
						!!! TODO: might neede to make this more flexible and make n_GLQ_cloudtop class ddependent for the differently sized classes
						! if (cloudtop_GLQ_point < n_GLQ_cloudtop) then
							! n1 = (GLQ_points_cloudtop(cloudtop_GLQ_point-1, class_number) + GLQ_points_cloudtop(cloudtop_GLQ_point, class_number)) / 2
							! n2 = (GLQ_points_cloudtop(cloudtop_GLQ_point, class_number) + GLQ_points_cloudtop(cloudtop_GLQ_point+1, class_number)) / 2
						! else
							!! print *, "this is not happening right?"
							! n1 = (GLQ_points_cloudtop(n_GLQ_cloudtop-1, class_number) + GLQ_points_cloudtop(n_GLQ_cloudtop, class_number)) / 2
							! n2 = n_clouds
						! end if
					! end if
					!!!
					n1 = cloudtop_GLQ_point
					n2 = cloudtop_GLQ_point
					
					! print *, "start through LWP indexes"
					! print *, "class_number", class_number
					! print *, "n1", n1
					! print *, "n2", n2
					do n = n1, n2
						! print *, "n", n
						fill_i = int(original_cloudtop_LWP_indexes(n, 1, class_number))
						fill_j = int(original_cloudtop_LWP_indexes(n, 2, class_number))
						
						! print *, "A"
						! print *, fill_i
						! print *, fill_j
						lwu(fill_i, fill_j,1:k1) =  lwUp_slice  (i,1:k1)
						lwd(fill_i, fill_j,1:k1) = -lwDown_slice(i,1:k1)
						! print *, "AA"
						!____________________
						if (.not. rad_longw) then !we get LW at surface identically to how it is done in sunray subroutine 
						!do i=2,i1
						lwd(fill_i, fill_j,1) =  -0.8 * boltz * thl0(fill_i, fill_j,1) ** 4.
						lwu(fill_i, fill_j,1) =  1.0 * boltz * tskin(fill_i, fill_j) ** 4.
						!end do
						end if
						!____________________
						! print *, "FFFFFF"
						swu(fill_i, fill_j,1:k1) =  swUp_slice  (i,1:k1)
						swd(fill_i, fill_j,1:k1) = -swDown_slice(i,1:k1)

						! print *, "BBBBBB"
						swdir(fill_i, fill_j,1:k1) = -swDownDir_slice(i,1:k1)
						swdif(fill_i, fill_j,1:k1) = -swDownDif_slice(i,1:k1)
						lwc  (fill_i, fill_j,1:k1) =  LWP_slice      (i,1:k1)

						! print *, "KKKKKK"
						lwuca(fill_i, fill_j,1:k1) =  lwUpCS_slice  (i,1:k1)
						lwdca(fill_i, fill_j,1:k1) = -lwDownCS_slice(i,1:k1)
						swuca(fill_i, fill_j,1:k1) =  swUpCS_slice  (i,1:k1)
						swdca(fill_i, fill_j,1:k1) = -swDownCS_slice(i,1:k1)

						! print *, "NNNNNN"
						SW_up_TOA (fill_i, fill_j) =  swUp_slice  (i,krad2)
						SW_dn_TOA (fill_i, fill_j) = -swDown_slice(i,krad2)
						LW_up_TOA (fill_i, fill_j) =  lwUp_slice  (i,krad2)
						LW_dn_TOA (fill_i, fill_j) = -lwDown_slice(i,krad2)

						! print *, "AAAAAA"
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