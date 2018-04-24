C
C  Copyright (C) 2009-2011, 2018  Smithsonian Astrophysical Observatory
C
C
C This program is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 3 of the License, or
C (at your option) any later version.
C
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
C

      subroutine check_bounds( n, par, lb, ub, epsfcn )
      implicit none
      integer n
      double precision par(n), lb(n), ub(n), epsfcn
      integer ii
      double precision eps, h, dpmpar
      
      eps = dsqrt( dmax1( epsfcn, dpmpar(1) ) )
      
      do ii = 1, n
         if ( par(ii) .lt. lb(ii) ) then
            h = eps * dabs( par( ii ) )
            if ( h .eq. 0.0d0 ) h = eps
            par(ii) = dmin1( lb(ii) + h, ub(ii) )
         end if
         if ( par(ii) .gt. ub(ii) ) then
            h = eps * dabs( par( ii ) )
            if ( h .eq. 0.0d0 ) h = eps
            par(ii) = dmax1( par(ii) - h, lb(ii) )
         end if
      enddo

      return
c
c     last card of subroutine check_bounds.
c
      end

      subroutine mylmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *     diag,mode,factor,nprint,info,nfev,fjac,ldfjac,ipvt,qtf,
     +     wa1,wa2,wa3,wa4,lb,ub,fmin)
c     *     lowtri,ifault,covarerr)
c     **********
c
c     subroutine mylmdif
      implicit none
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
      integer ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),
     *                 wa1(n),wa2(n),wa3(n),wa4(m)
      double precision lb(n),ub(n)
c      double precision lowtri(n*(n+1)/2)
      double precision fmin,enorm
      integer iflag, ii
      external fcn, covar
      iflag = 1

      if ( n .gt. m ) then
        print *, 'the number of parameters must < number of data points'
        info = 0
        return
      endif

      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *     diag,mode,factor,nprint,info,nfev,fjac,ldfjac,
     *     ipvt,qtf,wa1,wa2,wa3,wa4,lb,ub)
      fmin = enorm( m, fvec )**2.0
      call covar( n, fjac, ldfjac, ipvt, ftol, wa1 )
      return
c
c     last card of subroutine mylmdif.
c
      end
