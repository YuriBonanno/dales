module modradrrtmgyuri
implicit none

contains
	subroutine testyuri(tg_slice, cloudFrac, IWP_slice, LWP_slice, iceRe, liquidRe & !input
, tg_slice_reduced, cloudFrac_reduced, IWP_slice_reduced, LWP_slice_reduced, iceRe_reduced, liquidRe_reduced ) ! output
	
		use modglobal, only: imax,kind_rb

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
	end subroutine testyuri
	
end module modradrrtmgyuri