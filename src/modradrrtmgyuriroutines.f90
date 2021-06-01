module modradrrtmgyuriroutines

implicit none

contains

	!__________________________________________________________
	!__________________________________________________________
	!__________________________________________________________
	!__________________________________________________________

	! quicksort.f -*-f90-*-
	! Author: t-nissie
	! License: GPLv3
	! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
	!!
	recursive subroutine quicksortindexes(a, first, last, indexes)
	implicit none
	real*8  a(*), x, t
	integer first, last
	integer i, j
	integer indexes(last-first+1,2)
	integer temp_hor_index, temp_ver_index

	  x = a( (first+last) / 2 )
	  i = first
	  j = last
	  do
		 do while (a(i) < x)
			i=i+1
		 end do
		 do while (x < a(j))
			j=j-1
		 end do
		 if (i >= j) exit
		 t = a(i);  a(i) = a(j);  a(j) = t
		 temp_hor_index = indexes(i,1); indexes(i,1) = indexes(j,1); indexes(j,1) = temp_hor_index
		 temp_ver_index = indexes(i,2); indexes(i,2) = indexes(j,2); indexes(j,2) = temp_ver_index
		 i=i+1
		 j=j-1
	  end do
	  if (first < i-1) call quicksortindexes(a, first, i-1, indexes(first:i-1, :))
	  if (j+1 < last)  call quicksortindexes(a, j+1, last, indexes(j+1:last, :))
	end subroutine quicksortindexes

	! quicksort.f -*-f90-*-
	! Author: t-nissie
	! License: GPLv3
	! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
	!!
	recursive subroutine quicksort(a, first, last)
	implicit none
	real*8  a(*), x, t
	integer first, last
	integer i, j


	  x = a( (first+last) / 2 )
	  i = first
	  j = last
	  do
		 do while (a(i) < x)
			i=i+1
		 end do
		 do while (x < a(j))
			j=j-1
		 end do
		 if (i >= j) exit
		 t = a(i);  a(i) = a(j);  a(j) = t

		 i=i+1
		 j=j-1
	  end do
	  if (first < i-1) call quicksort(a, first, i-1)
	  if (j+1 < last)  call quicksort(a, j+1, last)
	end subroutine quicksort


	!+-------------------------------------------------------------------
		  SUBROUTINE gauleg(x1,x2,x,w,n)
	!+-------------------------------------------------------------------
		  INTEGER n
		  REAL x1,x2,x(n),w(n)
		  DOUBLE PRECISION EPS
		  PARAMETER (EPS=3.d-14)
		  INTEGER i,j,m

		  DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

		  m=(n+1)/2
		  xm=0.5d0*(x2+x1)
		  xl=0.5d0*(x2-x1)
		  do 12 i=1,m
			z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
	1       continue
			  p1=1.d0
			  p2=0.d0
			  do 11 j=1,n
				p3=p2
				p2=p1
				p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
	11        continue
			  pp=n*(z*p1-p2)/(z*z-1.d0)
			  z1=z

			  z=z1-p1/pp
			if(abs(z-z1).gt.EPS)goto 1
			x(i)=xm-xl*z
			x(n+1-i)=xm+xl*z
			w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
			w(n+1-i)=w(i)
	12    continue

		  return
		  END
		  
	!  +-----------------------------------------------------------------
	   subroutine quantiles (n_s, n_quantiles, std, x, q)

	   integer :: k, n_s, n_quantiles, index_11, index_12, index_q, index_median
	   real :: confidence_1
	   real :: x(n_s)
	   real :: q(n_quantiles)
	   logical :: std


	   confidence_1 = 0.684

	   index_11     = ((1.0 - confidence_1) / 2.0) * n_s
	   index_12     = n_s - index_11
	   index_median = n_s / 2.0

	   if (std) then

		  q(1) = x(index_11)
		  q(2) = x(index_median)
		  q(3) = x(index_12)

	   else

		  do k = 1, n_quantiles
			 index_q = (float(k) / float(n_quantiles + 1)) * n_s
			 q(k)    = x(index_q)
		  end do

	   end if

	   return
	   end

	subroutine writetofile(filename, dataset, dims)
	use modraddata
	use modglobal, only: imax, jmax, kmax, kind_rb

	implicit none

	!Files
		logical :: fileexists=.false.
		integer :: dims, i, j, k
		character(len = 16) :: filename
		real(kind=kind_rb) :: dataset (imax, jmax, kmax)

		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=filename, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			open(11, file=filename, status="old", position="append", action="write")
		else
			open(11, file=filename, status="new", action="write")
		end if
		
		if (dims == 1) then
			do i=1,imax
				write(11, *) dataset(i,1,1)
			end do
		end if
			
		if (dims == 2) then
			do i=1,imax
				do j=1,jmax
					write(11, *) dataset(i,j,1)
				end do
			end do
		end if
		
		if (dims == 3) then
			do i=1,imax
				do j=1,jmax
					do k=1,kmax
						write(11, *) dataset(i,j,k)
					end do
				end do
			end do
		end if

		close(11)


	end subroutine

end module modradrrtmgyuriroutines

