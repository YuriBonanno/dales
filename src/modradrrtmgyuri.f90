module modradrrtmgyuri
use modraddata

implicit none

contains
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

	subroutine testyuriLWP(LWP_collumns, LWP_flattened)
	
		use modglobal, only: imax, jmax, kmax, kind_rb, grav
		use modfields, only: ql0
		
		integer :: i,j,k
		logical :: fileexists=.false.
		character(len = 16) :: filename = 'LWPcontentfirsttest.txt'
		character(len = 40) :: writestring
	
		real(KIND=kind_rb) ::   LWP_collumns   (imax, jmax, krad1),       &
		                        layerMass      (imax, jmax, krad1),       &
								collumns_layerP   (imax, jmax, krad1),       &
								collumns_interfaceP   (imax, jmax, krad2),       &
								qcl_collumns   (imax, jmax, kradmax),     &
			                    LWP_flattened  (imax, jmax)
	
		qcl_collumns(:,:,:) = ql0(:,:,:)
	    do k=1, kmax
		   collumns_layerP(:,:,k) = presf_input(k)
		   collumns_interfaceP(:,:,p) = presh_input(k)
	    end do
	    collumns_interfaceP(:,:, krad2)  = min( 1.e-4_kind_rb , 0.25*collumns_layerP(1,:,krad1) )
	
	    do k=kmax+1, kradmax
		   collumns_layerP(:,:,k) = presf_input(k)
	    end do
	    collumns_layerP (:,:, krad1) = 0.5*presh_input(krad1)
	
	    do k=1, kradmax
		   collumns_layerMass(:,:,k) = 100.*( collumns_interfaceP(:,:,k) - collumns_interfaceP(:,:,k+1) ) / grav  !of full level
	       LWP_collumns(:,:,k) = qcl_collumns(:,:,k)*collumns_layerMass(:,:,k)*1e3
	    end do
	    collumns_layerMass(:,:,krad1) = 100.*( collumns_interfaceP(:,:,krad1) - collumns_interfaceP(:,:,krad2) ) / grav
        LWP_collumns(:,:,krad1) = 0.
	    !LWP_collumns is not cumulative I think, Also I think it has different values for different z
	    !I need to flatten it.
	    do i=1,imax
		   do j=1,jmax
		      LWP_flattened(i,j) = SUM(LWP_collumns(i,j,:))
		   end do
	    end do

		inquire(file=filename, exist=fileexists)
		print *, fileexists
		if (fileexists) then
			open(11, file=filename, status="old", position="append", action="write")
		else
			open(11, file=filename, status="new", action="write")
		end if
		do i=1,imax
		   do j=1,jmax
		      write(11, *) LWP_flattened(i,j)
		   end do
	    end do
		writestring = 'testyurifirst'
		write(11, *) writestring
		close(11)
	end subroutine testyuriLWP
	
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