module modradrrtmgyuri
use modraddata

implicit none

contains
!test Routine
	subroutine testlus

		logical :: fileexists=.false.
		character(len = 16) :: filename = 'testlus.txt'
		character(len = 40) :: writestring
	
		inquire(file=filename, exist=fileexists)
		print *, fileexists
		if (fileexists) then
			open(11, file=filename, status="old", position="append", action="write")
		else
			open(11, file=filename, status="new", action="write")
		end if
		writestring = 'ja testlus'
		write(11, *) writestring
		close(11)
	end subroutine testlus


!Important Routine
	subroutine testyuriLWP(LWP_collumns, LWP_flattened, cloudFrac)
	
		use modglobal, only: imax, jmax, kmax, kind_rb, grav
		use modfields, only: ql0
		
		integer :: i,j,k,inverse_k
		logical :: fileLWPexists=.false.
		logical :: fileql0exists=.false.
		character(len = 16) :: filenameLWP = 'LWPcontent.txt'
		character(len = 16) :: filenameql0 = 'ql0content.txt'
		character(len = 40) :: writestring
	
		real(KIND=kind_rb) ::   LWP_collumns   (imax, jmax, krad1),                &! LWP per gridpoint
		                        collumns_layerMass      (imax, jmax, krad1),       &! Mass per gridpoint
								collumns_layerP   (imax, jmax, krad1),             &! Pressure at the full gridpoint
								collumns_interfaceP   (imax, jmax, krad2),         &! Pressure at the half gridpoint
								qcl_collumns   (imax, jmax, kradmax),              &! cloud liquid content
			                    LWP_flattened  (imax, jmax),                       &! LWP flattened over Z direction
								cloudFrac      (imax, jmax),                       &! Fraction of clouds per collumn
								ztop_field     (imax, jmax)							! Height of the highest cloud
	
		!Define qcl
		qcl_collumns(:,:,:) = ql0(:,:,:)
		
		!Define Pressure
	    do k=1, kmax
		   collumns_layerP(:,:,k) = presf_input(k)
		   collumns_interfaceP(:,:,k) = presh_input(k)
	    end do
		do i=1,imax
	       collumns_interfaceP(i,:, krad2)  = min( 1.e-4_kind_rb , 0.25*collumns_layerP(1,:,krad1) )
		end do
	    
	    do k=kmax+1, kradmax
		   collumns_layerP(:,:,k) = presf_input(k)
	    end do
	    collumns_layerP (:,:, krad1) = 0.5*presh_input(krad1)
	
		!This one could be partially moved to the other k do loops
		!Define Layermass and LWP
	    do k=1, kradmax
		   collumns_layerMass(:,:,k) = 100.*( collumns_interfaceP(:,:,k) - collumns_interfaceP(:,:,k+1) ) / grav  !of full level
	       LWP_collumns(:,:,k) = qcl_collumns(:,:,k)*collumns_layerMass(:,:,k)*1e3
	    end do
	    collumns_layerMass(:,:,krad1) = 100.*( collumns_interfaceP(:,:,krad1) - collumns_interfaceP(:,:,krad2) ) / grav
        LWP_collumns(:,:,krad1) = 0.
		
		
	    !LWP_collumns is not cumulative I think, Also I think it has different values for different z
	    !I need to flatten it.
		!Define Ztop and cloudfrac
	    ztop_field(:,:) = 0.
		do i=1,imax
		   do j=1,jmax
			! Defines LWP_flattened and sets cloudFrac to 1 when LWP larger than zero
		      LWP_flattened(i,j) = SUM(LWP_collumns(i,j,:))
			  if (LWP_flattened(i,j)>0.0) then
				cloudFrac(i,j) = 1
			  end if
			  do k=1,kmax
			     inverse_k = kmax + 1 - k
				 !looks through the liquid in a gridpoint from top to bottom and assigns the first nonzero value to the ztop_field and then exits the vertical loop
				 if (ql0(i,j,inverse_k)>0.0) then
					ztop_field(i,j) = ql0(i,j,inverse_k) 
					EXIT
				 end if
			  end do
		   end do
	    end do

		!Select only the collumns with a nonzero cloudratio

		!Order on basis of ztop_field
		
		!Order on basis of LWP



		inquire(file=filenameLWP, exist=fileLWPexists)
		inquire(file=filenameql0, exist=fileql0exists)
		!print *, fileLWPexists
		!print *, fileql0exists
		if (fileLWPexists) then
			open(11, file=filenameLWP, status="old", position="append", action="write")
		else
			open(11, file=filenameLWP, status="new", action="write")
		end if
		do i=1,imax
		   do j=1,jmax
		      write(11, *) LWP_flattened(i,j)
		   end do
	    end do
		writestring = 'LWPOneSet'
		write(11, *) writestring
		close(11)
		
		if (fileql0exists) then
			open(12, file=filenameql0, status="old", position="append", action="write")
		else
			open(12, file=filenameql0, status="new", action="write")
		end if
		do i=1,imax
		   do j=1,jmax
			  do k=1,kmax
				write(12, *) ql0(i,j,k)
			  end do
		   end do
	    end do
		writestring = 'ql0OneSet'
		write(12, *) writestring
		close(12)
		
	end subroutine testyuriLWP
	
	!OutdatedRoutine
	subroutine testyurirad(tg_slice, cloudFrac, IWP_slice, LWP_slice, iceRe, liquidRe & !input
, tg_slice_reduced, cloudFrac_reduced, IWP_slice_reduced, LWP_slice_reduced, iceRe_reduced, liquidRe_reduced ) ! output
	
		use modglobal, only: imax, kind_rb

		logical :: fileexists=.false.
		character(len = 16) :: filename = 'slicecontent.txt'
		character(len = 40) :: writestring
	

		!krad1 is a SAVE'd value in modraddata
		!real(KIND=kind_rb),intent(out) ::    LWP_slice(imax,krad1), &
		!                                     IWP_slice(imax,krad1), &
		!                                     cloudFrac(imax,krad1), &
		!                                     liquidRe (imax,krad1), &
		!                                     iceRe    (imax,krad1), &
		!									 tg_slice (imax)
	
		real(KIND=kind_rb) ::    LWP_slice(imax,krad1), &
								IWP_slice(imax,krad1), &
								cloudFrac(imax,krad1), &
								liquidRe (imax,krad1), &
								iceRe    (imax,krad1), &
								tg_slice (imax)
										 
		real(KIND=kind_rb) ::    LWP_slice_reduced(imax,krad1), &
								IWP_slice_reduced(imax,krad1), &
								cloudFrac_reduced(imax,krad1), &
								liquidRe_reduced (imax,krad1), &
								iceRe_reduced    (imax,krad1), &
								tg_slice_reduced (imax)
	
	
		inquire(file=filename, exist=fileexists)
		print *, fileexists
		if (fileexists) then
			open(11, file=filename, status="old", position="append", action="write")
		else
			open(11, file=filename, status="new", action="write")
		end if
		writestring = 'abcdefghijklmnopqrstuvwxyz'
		write(11, *) writestring
		close(11)
	end subroutine testyurirad
	
end module modradrrtmgyuri