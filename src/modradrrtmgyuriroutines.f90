module modradrrtmgyuriroutines
use modglobal, only: kind_rb, runtime
use modraddata, only: tnext
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
	recursive subroutine quicksortindexes(a, first, last, indexes, length)
	implicit none
	integer first, last, length
	integer i, j
	integer temp_hor_index, temp_ver_index
	integer indexes(length, 2)
	real(kind=kind_rb), dimension(:) :: a(length)
	real(kind=kind_rb) :: x, t
	
	!Determine the middle of the class
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
		 !Also move the indexes
		 temp_hor_index = indexes(i,1); indexes(i,1) = indexes(j,1); indexes(j,1) = temp_hor_index
		 temp_ver_index = indexes(i,2); indexes(i,2) = indexes(j,2); indexes(j,2) = temp_ver_index
		 i=i+1
		 j=j-1
	  end do

	!Recursively go through this subroutine to order the array
	  if (first < i-1) call quicksortindexes(a, first, i-1, indexes, length)
	  if (j+1 < last)  call quicksortindexes(a, j+1, last, indexes, length)
	end subroutine quicksortindexes


!This function determines the Gauss-Legendre Quadrature points
!+-------------------------------------------------------------------
  subroutine gauleg(x1,x2,x,w,n)
!+-------------------------------------------------------------------
    implicit none
    INTEGER n
    REAL(kind=kind_rb) :: x1,x2,x(n),w(n)
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
  end subroutine gauleg
		  
		  
	! This function determines the quantiles for the clouded classes, these quantiles are evenly spaced over the amount of points
	!  +-----------------------------------------------------------------
  subroutine quantiles (n_s, n_quantiles, std, x, q)
   implicit none
   integer :: k, n_s, n_quantiles, index_11, index_12, index_q, index_median
   real(kind=kind_rb) :: confidence_1
   real(kind=kind_rb) :: x(n_s)
   real(kind=kind_rb) :: q(n_quantiles)
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
     index_q = nint((float(k) / float(n_quantiles + 1)) * n_s)
     q(k)    = x(index_q)

    end do

   end if

   return
  end subroutine quantiles


	subroutine GetDiff(dataset, resultDataSet, xsize, ysize, zsize)
		integer :: xsize, ysize, zsize
		integer :: dims, i, j, k, m
		real(kind=kind_rb) :: dataset (xsize, ysize, zsize)
		real(kind=kind_rb) :: resultDataSet (xsize, ysize, zsize-1)
	
		do i=1,xsize
			do j=1,ysize
				do k=1,zsize-1
					resultDataSet(i,j,k) = dataset(i,j,k+1) - dataset(i,j,k)
				enddo
			enddo
		enddo
	
	end subroutine GetDiff
	! --------------------------------------------------------------------
	! PROGRAM  MeanVariance:
	!    This program reads in an unknown number of real values and
	! computes its mean, variance and standard deviation.  It contains
	! three subroutines:
	!    (1)  Sums()     - computes the sum and sum of squares of the input
	!    (2)  Result()   - computes the mean, variance and standard
	!                      deviation from the sum and sum of squares
	!    (3)  PrintResult() - print results
	! --------------------------------------------------------------------
	subroutine MeanVariance(dataset, filename, xsize, ysize, zsize)
		integer :: xsize, ysize, zsize
		REAL(kind=kind_rb)    :: SumMean, SumVar, x
		REAL(kind=kind_rb)    :: Mean(zsize), Std(zsize), Var(zsize)
		integer :: Ncolumns
		integer :: dims, i, j, k, m
		real(kind=kind_rb) :: dataset (xsize, ysize, zsize)
		logical :: fileexists=.false.
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		character(5000) :: frmt

		makedir = "datadir"
		call execute_command_line ('mkdir -p ' // trim(makedir))
		! for mean
		fullpath = trim(makedir) // '/' // trim(filename) // '_mean'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(10, file=fullpath, status="old", position="append", action="write")
		else
			open(10, file=fullpath, status="new", action="write")
		end if

		! for stdev
		fullpath = trim(makedir) // '/' // trim(filename) // '_stdev'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(11, file=fullpath, status="old", position="append", action="write")
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		! for var
		fullpath = trim(makedir) // '/' // trim(filename) // '_var'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(12, file=fullpath, status="old", position="append", action="write")
		else
			open(12, file=fullpath, status="new", action="write")
		end if

		Ncolumns = xsize*ysize
		x = 0.0
		SumMean = 0.0
		SumVar = 0.0

		frmt = "(F18.10"
		do k=1,zsize
			x = 0.0
			SumMean = 0.0
			SumVar = 0.0
			do i=1,xsize
				do j=1,ysize
					x = dataset(i, j, k)
					SumMean = SumMean + x
					SumVar = SumVar + x*x
				enddo
			enddo
			call  Results(SumMean, SumVar, Ncolumns, Mean(k), Std(k), Var(k))  ! compute results
			if (k==1) then
				continue
			else
				frmt = trim(frmt)
				frmt = trim(frmt) // ",F18.10 "
				frmt = trim(frmt)
			end if

		enddo
		frmt = trim(frmt) // ", A)"
		
		write(10, frmt, advance="no")  Mean(:), ","
		write(11, frmt, advance="no")  Std(:), ","
		write(12, frmt, advance="no")  Var(:), ","
		close(10)
		close(11)
		close(12)
		! deallocate(fullpath)
		! deallocate(makedir)
	end subroutine MeanVariance

	subroutine MeanVarianceOnlyClouds(dataset, filename, xsize, ysize, zsize, NColumns)
		integer :: xsize, ysize, zsize
		REAL(kind=kind_rb)    :: SumMean, SumVar, x
		REAL(kind=kind_rb)    :: Mean(zsize), Std(zsize), Var(zsize)
		integer :: Ncolumns
		integer :: dims, i, j, k, m
		real(kind=kind_rb) :: dataset (xsize, ysize, zsize)
		logical :: fileexists=.false.
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		character(5000) :: frmt

		makedir = "datadir"
		call execute_command_line ('mkdir -p ' // trim(makedir))
		! for mean
		fullpath = trim(makedir) // '/' // trim(filename) // '_mean'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(10, file=fullpath, status="old", position="append", action="write")
		else
			open(10, file=fullpath, status="new", action="write")
		end if

		! for stdev
		fullpath = trim(makedir) // '/' // trim(filename) // '_stdev'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(11, file=fullpath, status="old", position="append", action="write")
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		! for var
		fullpath = trim(makedir) // '/' // trim(filename) // '_var'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(12, file=fullpath, status="old", position="append", action="write")
		else
			open(12, file=fullpath, status="new", action="write")
		end if

		x = 0.0
		SumMean = 0.0
		SumVar = 0.0

		if (NColumns /= 0) then
			frmt = "(F18.10"
			do k=1,zsize
				x = 0.0
				SumMean = 0.0
				SumVar = 0.0
				do i=1,xsize
					do j=1,ysize
						x = dataset(i, j, k)
						SumMean = SumMean + x
						SumVar = SumVar + x*x
					enddo
				enddo
				call  Results(SumMean, SumVar, Ncolumns, Mean(k), Std(k), Var(k))  ! compute results
				if (k==1) then
					continue
				else
					frmt = trim(frmt)
					frmt = trim(frmt) // ",F18.10 "
					frmt = trim(frmt)
				end if

			enddo
			frmt = trim(frmt) // ", A)"
		else
			Mean(:) = 0
			Std(:) = 0
			Var(:) = 0
			frmt = "(F18.10"
			do k=1,zsize
				if (k==1) then
					continue
				else
					frmt = trim(frmt)
					frmt = trim(frmt) // ",F18.10 "
					frmt = trim(frmt)
				end if
			end do
			frmt = trim(frmt) // ", A)"
		end if
		
		write(10, frmt, advance="no")  Mean(:), ","
		write(11, frmt, advance="no")  Std(:), ","
		write(12, frmt, advance="no")  Var(:), ","
		close(10)
		close(11)
		close(12)
		! deallocate(fullpath)
		! deallocate(makedir)
	end subroutine MeanVarianceOnlyClouds
	
	
	! --------------------------------------------------------------------
	! SUBROUTINE  Results():
	!    This subroutine computes the mean, variance and standard deviation
	! from the sum and sum-of-squares:
	!    (1) Sum       - sum of input values
	!    (2) SumSQR    - sun-of-squares
	!    (3) n         - number of input data items
	!    (4) Mean      - computed mean value
	!    (5) Variance  - computed variance
	!    (6) StdDev    - computed standard deviation
	! --------------------------------------------------------------------
	SUBROUTINE  Results(Sum1, SumSQR, n, Mean, StdDev, Variance)
		IMPLICIT  NONE

		INTEGER, INTENT(IN) :: n
		REAL(kind=kind_rb), INTENT(IN)    :: Sum1, SumSQR
		REAL(kind=kind_rb), INTENT(OUT)   :: Mean, StdDev, Variance

		Mean = Sum1 / n
		Variance = (SumSQR - Sum1*Sum1/n)/(n-1)
		if (Variance < 0.0) then
			Variance = 0
		end if
		StdDev   = SQRT(Variance)
	END SUBROUTINE

	subroutine finishstatisticsline(filename)
		logical :: fileexists=.false.
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
	
	
		makedir = "datadir"
		call execute_command_line ('mkdir -p ' // trim(makedir))
		! for mean
		fullpath = trim(makedir) // '/' // trim(filename) // '_mean'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(10, file=fullpath, status="old", position="append", action="write")
		else
			open(10, file=fullpath, status="new", action="write")
		end if

		! for stdev
		fullpath = trim(makedir) // '/' // trim(filename) // '_stdev'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(11, file=fullpath, status="old", position="append", action="write")
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		! for var
		fullpath = trim(makedir) // '/' // trim(filename) // '_var'
		inquire(file=fullpath, exist=fileexists)
		if (fileexists) then
			open(12, file=fullpath, status="old", position="append", action="write")
		else
			open(12, file=fullpath, status="new", action="write")
		end if
	
	
		write(10, *)  ' '
		write(11, *)  ' '
		write(12, *)  ' '

		close(10)
		close(11)
		close(12)
		! deallocate(fullpath)
		! deallocate(makedir)
	
	end subroutine finishstatisticsline
	
	! write to file with real(kind=kind_rb) with assumed size imax, jmax, kradmax
	subroutine writetofile(filename, dataset, dims, printLast)
	use modraddata
	use modglobal, only: imax, jmax, kmax, kind_rb

	! implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: dims, i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		real(kind=kind_rb) :: dataset (imax, jmax, kradmax)
		character(50000) :: frmt

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if
		
		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		frmt = "(F18.10"
		do i=2,imax
			frmt = trim(frmt)
			frmt = trim(frmt) // ",F18.10 "
			frmt = trim(frmt)
		end do
		frmt = trim(frmt) // ")"
		
		if (dims == 1) then
			!do i=1,imax
				write(11, frmt) dataset(:,1,1)
			!end do
		end if
			
		if (dims == 2) then
			!do i=1,imax
				do j=1,jmax
					write(11, frmt) dataset(:,j,1)
				end do
			!end do
		end if
		
		if (dims == 3) then
			!do i=1,imax
				do j=1,jmax
					do k=1,kmax
						write(11, frmt) dataset(:,j,k)
					end do
				end do
			!end do
		end if

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine

	! write to file with real(kind=kind_rb) with defined size
	subroutine writetofiledefinedsize(filename, dataset, dims, xsize, ysize, zsize, printLast)
	! use modraddata
	use modglobal, only: kind_rb

	! implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: dims, i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		integer :: xsize, ysize, zsize
		real(kind=kind_rb) :: dataset (xsize, ysize, zsize)
		character(50000) :: frmt

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		! print *, "test start"
		
		! print *, "now starting the concatenation"
		frmt = "(F18.10"
		do i=2,xsize
			frmt = trim(frmt)
			frmt = trim(frmt) // ",F18.10 "
			frmt = trim(frmt)
		end do
		frmt = trim(frmt) // ")"
		
		! print *, "finished concatenation"
		
		! print *, "concatenation result:"
		! print *, frmt
		! print *, "trimmed concatenation result:"
		! print *, trim(frmt)
		
		! print *, "now writing to datadir"
		if (dims == 1) then
			!print *, "1D write"
			! do i=1,imax
				write(11, frmt) dataset(:,1,1)
			! end do
		end if
			
		if (dims == 2) then
			! do i=1,imax
				do j=1,ysize
					write(11, frmt) dataset(:,j,1)
				end do
			! end do
		end if
		
		if (dims == 3) then
			! do i=1,imax
				do j=1,ysize
					do k=1,zsize
						write(11, frmt) dataset(:,j,k)
					end do
				end do
			! end do
		end if

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine
	
		! write to file with real(kind=kind_rb) with defined size
	subroutine writetofiledefinedsizedifferentlengths(filename, dataset, sizeset, dims, xsize, ysize, zsize, printLast)
	! use modraddata
	use modglobal, only: kind_rb

	! implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: dims, i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		integer :: xsize, ysize, zsize
		integer, dimension(:) :: sizeset
		real(kind=kind_rb) :: dataset (xsize, ysize, zsize)
		character(50000) :: frmt

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		! print *, "test start"
		
		! print *, "now starting the concatenation"
		frmt = "(F18.10"
		do i=2,xsize
			frmt = trim(frmt)
			frmt = trim(frmt) // ",F18.10 "
			frmt = trim(frmt)
		end do
		frmt = trim(frmt) // ")"
		
		! print *, "finished concatenation"
		
		! print *, "concatenation result:"
		! print *, frmt
		! print *, "trimmed concatenation result:"
		! print *, trim(frmt)
		
		! print *, "now writing to datadir"
		if (dims == 1) then
			!print *, "1D write"
			! do i=1,imax
				write(11, frmt) dataset(:,1,1)
			! end do
		end if
			
		if (dims == 2) then
			! do i=1,imax
				do j=1,ysize
					write(11, frmt) dataset(:,j,1)
				end do
			! end do
		end if
		
		if (dims == 3) then
			! do i=1,imax
				do j=1,ysize
					do k=1,zsize
						write(11, frmt) dataset(:,j,k)
					end do
				end do
			! end do
		end if

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)
	end subroutine
	
	! write to file with integers with defined size
	subroutine writetofiledefinedsizeint(filename, dataset, dims, xsize, ysize, zsize, printLast)
	! use modraddata
	! use modglobal, only: imax, jmax, kmax, kind_rb

	implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: dims, i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		integer :: xsize, ysize, zsize
		integer :: dataset (xsize, ysize, zsize)
		character(50000) :: frmt

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		frmt = "(I6"
		do i=2,xsize
			frmt = trim(frmt)
			frmt = trim(frmt) // ",I6 "
			frmt = trim(frmt)
		end do
		frmt = trim(frmt) // ")"
		
		if (dims == 1) then
			!do i=1,imax
				write(11, frmt) dataset(:,1,1)
			!end do
		end if
			
		if (dims == 2) then
			!do i=1,imax
				do j=1,ysize
					write(11, frmt) dataset(:,j,1)
				end do
			!end do
		end if
		
		if (dims == 3) then
			!do i=1,imax
				do j=1,ysize
					do k=1,zsize
						write(11, frmt) dataset(:,j,k)
					end do
				end do
			!end do
		end if

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine

	! write single ints to file
	subroutine writeinttofile(filename, intvalue, printLast)
	! use modraddata

	implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		integer :: intvalue

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
	
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		write(11, *) intvalue

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine

	! write single real to file
	subroutine writerealtofile(filename, realvalue, printLast)
	! use modraddata

	implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		real(kind=kind_rb) :: realvalue

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
	
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		write(11, *) realvalue

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine

	! subroutine handleIndexes(dataset, filename, length)
		! integer :: length
		! integer :: dims, i
		! real(kind=kind_rb) :: dataset (length)
		! logical :: fileexists=.false.
		! character(*) :: filename
		! character(:), allocatable :: fullpath
		! character(:), allocatable :: makedir
		! character(500) :: frmt

		! makedir = "datadir"
		! call execute_command_line ('mkdir -p ' // trim(makedir))
		! ! ! ! for mean
		! fullpath = trim(makedir) // '/' // trim(filename)'
		! inquire(file=fullpath, exist=fileexists)
		! if (fileexists) then
			! open(10, file=fullpath, status="old", position="append", action="write")
		! else
			! open(10, file=fullpath, status="new", action="write")
		! end if

		! frmt = "("
		! do i=1,length
			! frmt = trim(frmt)
			! frmt = trim(frmt) // ",F18.10 "
			! frmt = trim(frmt)
			
		! enddo
		! frmt = trim(frmt) // ", A)"#
		
		
		! write to file with real(kind=kind_rb) with defined size
	subroutine testwritetofiledefinedsize(filename, dataset, dims, xsize, ysize, zsize, printLast)
	! use modraddata
	use modglobal, only: kind_rb

	! implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: dims, i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		integer :: xsize, ysize, zsize
		real(kind=kind_rb) :: dataset (xsize, ysize, zsize)
		character(50000) :: frmt

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		! print *, "test start"
		
		! print *, "now starting the concatenation"
		frmt = "(F18.10"
		do i=2,xsize
			frmt = trim(frmt)
			frmt = trim(frmt) // ",F18.10 "
			frmt = trim(frmt)
		end do
		frmt = trim(frmt) // ")"
		
		! print *, "finished concatenation"
		
		! print *, "concatenation result:"
		! print *, frmt
		! print *, "trimmed concatenation result:"
		! print *, trim(frmt)
		
		! print *, "now writing to datadir"
		if (dims == 1) then
			!print *, "1D write"
			! do i=1,xsize
				write(11, frmt) dataset(:,1,1)
			! end do
		end if
			
		if (dims == 2) then
			do j=1,ysize
				! do i=1,xsize
					write(11, frmt) dataset(:,j,1)
				! end do
			end do
		end if
		
		if (dims == 3) then
			do k=1,zsize
				do j=1,ysize
					! do i=1,xsize
						write(11, frmt) dataset(:,j,k)
					! end do
				end do
			end do
		end if

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine
		
		! write(10, frmt, advance="no")  dataset(:), ","
		! close(10)
		! deallocate(fullpath)
		! deallocate(makedir)
	! end subroutine MeanVariance

	subroutine testwritetofiledefinedsizeint(filename, dataset, dims, xsize, ysize, zsize, printLast)
	! use modraddata
	! use modglobal, only: imax, jmax, kmax, kind_rb

	implicit none

	!Files
		logical :: fileexists=.false.
		logical :: printLast
		integer :: dims, i, j, k
		character(*) :: filename
		character(200) :: fullpath
		character(200) :: makedir
		integer :: xsize, ysize, zsize
		integer :: dataset (xsize, ysize, zsize)
		character(50000) :: frmt

		if (printLast) then
			if (.NOT. ((tnext/1000)>=runtime)) then
				return
			end if
		end if

		makedir = "datadir"
		!call execute_command_line ('mkdir -p out/' // adjustl(trim( makedir ) ) )
		call execute_command_line ('mkdir -p ' // trim(makedir))
		fullpath = trim(makedir) // '/' // trim(filename)
		
		!__________________________________________________________
		!Make and write to files
		if (dims<1 .or. dims>3) then
			print *, "writing to file failed because of erroneous dimensions"
			return
		end if
		
		inquire(file=fullpath, exist=fileexists)
		!print *, fileexists
		if (fileexists) then
			! open(11, file=fullpath, status="old", position="append", action="write")
			if (printLast) then
				open(11, file=fullpath, status="replace", action="write")
			else
				open(11, file=fullpath, status="old", position="append", action="write")
			end if
		else
			open(11, file=fullpath, status="new", action="write")
		end if
		
		frmt = "(I6"
		do i=2,xsize
			frmt = trim(frmt)
			frmt = trim(frmt) // ",I6 "
			frmt = trim(frmt)
		end do
		frmt = trim(frmt) // ")"
		
		if (dims == 1) then
			!do i=1,imax
				write(11, frmt) dataset(:,1,1)
			!end do
		end if
			
		if (dims == 2) then
			do j=1,ysize
				! do i=1,xsize
					write(11, frmt) dataset(:,j,1)
				! end do
			end do
		end if
		
		if (dims == 3) then
			do k=1,zsize
				do j=1,ysize
					! do i=1,xsize
						write(11, frmt) dataset(:,j,k)
					! end do
				end do
			end do
		end if

		close(11)
		! deallocate(fullpath)
		! deallocate(makedir)

	end subroutine

end module modradrrtmgyuriroutines


