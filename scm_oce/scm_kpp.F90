module scm_kpp
!
   implicit none
!
contains


!===================================================================================================
      SUBROUTINE lmd_vmix(N,Akv,Akt,u,v,rho1,bvf,z_r)
!---------------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in   )                   :: N
      real(8), intent(  out)                   :: Akv (0:N  )
      real(8), intent(  out)                   :: Akt (0:N,2)
      real(8), intent(in   )                   :: u   (1:N  )
      real(8), intent(in   )                   :: v   (1:N  )
      real(8), intent(in   )                   :: rho1(1:N  )
      real(8), intent(in   )                   :: z_r (1:N  )
      real(8), intent(in   )                   :: bvf (0:N  )
! local variables
      integer                                  :: k
      real(8)                                  :: Rig(1:N)
      real(8)                                  :: nu_sx,dudz,dvdz,cff
      real(8), parameter                       :: eps = 1.E-16
      real(8), parameter :: nuwm = 1.0e-4     !<-- minimum turbulent viscosity
      real(8), parameter :: nuws = 0.1e-4     !<-- minimum turbulent diffusion
      real(8), parameter :: nu0m = 50.e-4     !<-- KPP parameter
      real(8), parameter :: nu0s = 50.e-4     !<-- KPP parameter
      real(8), parameter :: nu0c = 0.1        !<-- turbulent viscosity/diffusion when convective adjustement
      real(8), parameter :: Ri0  = 0.7
!---------------------------------------------------------------------------------------------------
! Compute interior viscosities and diffusivities everywhere
! as the superposition of:
! (i) internal wave breaking
! (ii) local Richardson number instability due to resolved
!      vertical shear and convection due to static instability
!-------
!
! Richardson number mixing & convective adjustment
!-------
      do k=1,N-1
        cff    =  1./ ( z_r(k+1) - z_r(k) )
        dudz   =  cff*(   u(k+1) -   u(k) )       ! vertical shear
        dvdz   =  cff*(   v(k+1) -   v(k) )
        Rig(k) = bvf(k)/(Ri0*max(dudz**2+dvdz**2,eps)) ! Gradient Richardson number
      enddo

      do k=1,N-1
        cff   = min( 1., max(0.,Rig(k)) )
        nu_sx = 1.-cff*cff
        nu_sx = nu_sx*nu_sx*nu_sx           ! mixing due to vertical
        AKv(k  )=nuwm +nu0m*nu_sx           ! shear instability and
        AKt(k,1)=nuws +nu0s*nu_sx           ! internal wave breaking

        if (rho1(k+1).gt.rho1(k)) then          ! interior convective
          AKv(k  )=AKv(k      )+nu0c        ! diffusivity due to
          AKt(k,1)=AKt(k,1)+nu0c        ! static instability
        endif
        Akt(k,2)=AKt(k,1)
      enddo                   ! <-- k

      AKv(0  )=0.
      AKv(N  )=0.
      AKt(0,1)=0.
      AKt(N,1)=0.
      AKt(0,2)=0.
      AKt(N,2)=0.
!---------------------------------------------------------------------------------------------------
end subroutine lmd_vmix
!===================================================================================================




!===================================================================================================
      SUBROUTINE lmd_kpp(N,Akv,Akt,hbl,u,v,t,bvf,z_r,z_w,Hz,Ricr,f,sustr,svstr,srflx,stflx,rho0,swr_frac,ghat)
!---------------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in   )                   :: N
      real(8), intent(  out)                   :: Akv     (0:N  )
      real(8), intent(  out)                   :: Akt     (0:N,2)
      real(8), intent(in   )                   :: u       (1:N  )
      real(8), intent(in   )                   :: v       (1:N  )
      real(8), intent(in   )                   :: t       (1:N,2)
      real(8), intent(in   )                   :: bvf     (0:N  )
      real(8), intent(in   )                   :: z_r     (1:N  )
      real(8), intent(in   )                   :: Hz      (1:N  )
      real(8), intent(in   )                   :: z_w     (0:N  )
      real(8), intent(in   )                   :: swr_frac(0:N  )
      real(8), intent(inout)                   :: ghat    (0:N  )
      real(8), intent(inout)                   :: hbl
      real(8), intent(in   )                   :: Ricr,f,rho0
      real(8), intent(in   )                   :: sustr,svstr
      real(8), intent(in   )                   :: srflx,stflx(2)
! local variables
      integer                                  :: k,kbl
      real(8)     :: Cg,Vtc,Cr(0:N),FC(0:N),z_bl,zscale
      real(8)     :: du(0:N),dv(0:N),Bfsfc,Bo,Bosol,Vtsq
      real(8)     :: ws,wm,cff_up,cff_dn,f1,sigma,Kern
      real(8)     :: a1,a2,a3,alpha,beta,Ri_inv,cff,ustar
      real(8)     :: Av_bl,dAv_bl,Gm1,dGm1dS
      real(8)     :: At_bl,dAt_bl,Gt1,dGt1dS
      real(8)     :: As_bl,dAs_bl,Gs1,dGs1dS
      real(8)     :: bvf_r(1:N)
      real(8), parameter :: epssfc =  0.1
      real(8), parameter :: betaT  = -0.2
      real(8), parameter :: c_s    = 98.96
      real(8), parameter :: C_Ek   = 258.
      real(8), parameter :: Cstar  = 10.
      real(8), parameter :: Cv     = 1.8
      real(8), parameter :: vonKar = 0.4
      real(8), parameter :: eps    = 1.D-14
      real(8), parameter :: g      =  9.81
!---------------------------------------------------------------------------------------------------
!
       Ri_inv= 1./Ricr
!
! Nondimensional constants for computing non-local flux and convective
! deepening of surface boundary layer.
!------
      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)
!
! Compute thermal expansion, "alpha" [kg/m^3/degC], and saline
! contraction, "beta" [kg/m^3/PSU], coefficients at surface
!------
      call alfabeta(N,alpha,beta,t,rho0)
!
! compute surface turbulent buoyancy forcing "Bo" [m^2/s^3] (in doing
! so remove incoming solar shortwave radiation component and save it
! separately as "Bosol")
!------
      Bo   =g*(alpha*(stflx(1)-srflx)-beta * stflx(2))
      Bosol=g* alpha* srflx
!
! compute turbulent  friction velocity "ustar" from wind stress
!------
      ustar=sqrt( sqrt(sustr**2+svstr**2) )

!      hbl = hbl_old           !<-- boundary layer depth from previous time step
      kbl = 0
      Cr(0)=0.
      Cr(N)=0.
      FC(N)=0.

!
! compute vertical shear
!------
      do k=1,N-1
          cff    =  1./(z_r(k+1)-z_r(k))
          du(k)  = cff*(  u(k+1)-  u(k))
          dv(k)  = cff*(  v(k+1)-  v(k))
      enddo

      du(N)=du(N-1)
      dv(N)=dv(N-1)
      du(0)=du(  1)
      dv(0)=dv(  1)

!
! Integral condition to find hbl
!------
      do k=N,1,-1

        zscale = z_w(N)-z_w(k-1)
        Kern   = zscale/(zscale+epssfc*hbl)
        Bfsfc  = Bo +Bosol*(1.-swr_frac(k-1))
        Call lmd_wscale(Bfsfc,zscale,ustar,wm,ws,hbl)

        cff    = bvf(k)*bvf(k-1)
        if (cff.gt.0.D0) then
            cff=cff/(bvf(k)+bvf(k-1))
        else
            cff=0.D0
        endif

        bvf_r(k) = cff + 0.25*(bvf(k)+bvf(k-1))

        FC(k-1)=FC(k) + Kern*Hz(k)*( 0.375*(                         &
               du(k)**2+du(k-1)**2+dv(k)**2+dv(k-1)**2)              &
                   +0.25 *(du(k-1)*du(k)+dv(k-1)*dv(k))              &
                    -Ri_inv*( cff + 0.25*(bvf(k)+bvf(k-1)))          &
                                           - C_Ek*f*f      )

        Vtsq=Vtc*ws*sqrt(max(0., bvf(k-1)))
        Cr(k-1)=FC(k-1) +Vtsq
        if (kbl.eq.0 .and.  Cr(k-1).lt.0.) kbl=k

      enddo

!
! From kbl determine hbl by interpolation
!-------
      if (kbl.gt.0) then
          k=kbl
          hbl=z_w(N)- ( z_w(k-1)*Cr(k)-z_w(k)*Cr(k-1)    &
                                      )/(Cr(k)-Cr(k-1))
      else
          hbl=z_w(N)-z_w(0)+eps
      endif

!
! Find buoyancy forcing for final "hbl" values
!-------
      k=kbl
      z_bl=z_w(N)-hbl
      zscale=hbl

      if (k.gt.0 .and. swr_frac(k-1).gt. 0.) then
          Bfsfc=Bo +Bosol*( 1. -swr_frac(k-1)         &
                  *swr_frac(k)*(z_w(k)-z_w(k-1))      &
                   /( swr_frac(k  )*(z_w(k)   -z_bl)  &
                     +swr_frac(k-1)*(z_bl -z_w(k-1))  &
                                                           ))
      else
          Bfsfc=Bo+Bosol
      endif
!
! Compute tubulent velocity scales (wm,ws) at "hbl"
!-------
       Call lmd_wscale(Bfsfc,zscale,ustar,wm,ws,hbl)
!
! Compute nondimensional shape function coefficients Gx( ) by
! matching values and vertical derivatives of interior mixing
! coefficients at hbl (sigma=1).
!-------


      f1=5.0 * max(0., Bfsfc) * vonKar/(ustar**4+eps)

      cff=1./(z_w(k)-z_w(k-1))
      cff_up=cff*(z_bl -z_w(k-1))
      cff_dn=cff*(z_w(k)  -z_bl)

      Av_bl=cff_up*AKv(k)+cff_dn*AKv(k-1)
      dAv_bl=cff * (AKv(k)   -   AKv(k-1))
      Gm1=Av_bl/(hbl*wm+eps)
      dGm1dS=min(0., Av_bl*f1-dAv_bl/(wm+eps))

      At_bl =cff_up*AKt(k,1)+cff_dn*AKt(k-1,1)
      dAt_bl=cff * (AKt(k,1)-AKt(k-1,1))
      Gt1=At_bl/(hbl*ws+eps)
      dGt1dS=min(0., At_bl*f1-dAt_bl/(ws+eps))

      As_bl =cff_up*AKt(k,2)+cff_dn*AKt(k-1,2)
      dAs_bl=cff * (AKt(k,2)-AKt(k-1,2))
      Gs1=At_bl/(hbl*ws+eps)
      dGs1dS=min(0., As_bl*f1-dAs_bl/(ws+eps))

! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.

       ghat(0:N) = 0.

       do k=N-1,kbl,-1

         zscale=z_w(N)-z_w(k)
         Call lmd_wscale(Bfsfc,zscale,ustar,wm,ws,hbl)

! Compute vertical mixing coefficients
         sigma=(z_w(N)-z_w(k))/max(hbl,eps)

         a1=sigma-2.
         a2=3.-2.*sigma
         a3=sigma-1.

         if (sigma.lt.0.07D0) then
            cff=0.5*(sigma-0.07D0)**2/0.07D0
         else
            cff=0.D0
         endif

         AKv(k      )=wm*hbl*( cff + sigma*( 1.+sigma*(    &
                                 a1+a2*Gm1+a3*dGm1dS )))
         AKt(k,1)=ws*hbl*( cff + sigma*( 1.+sigma*(    &
                                 a1+a2*Gt1+a3*dGt1dS )))
         AKt(k,2)=ws*hbl*( cff + sigma*( 1.+sigma*(    &
                                 a1+a2*Gs1+a3*dGs1dS )))
! Compute non-local term
         if (Bfsfc .lt. 0.) then
            ghat(k)=Cg * sigma*(1.-sigma)**2
            !ghat(k)=Cg * sigma*( 1.+sigma*( a1+a2*Gm1+a3*dGm1dS ))
         else
            ghat(k)=0.
         endif

       enddo

!---------------------------------------------------------------------------------------------------
END SUBROUTINE lmd_kpp
!===================================================================================================
!








!===================================================================================================
      subroutine lmd_wscale (Bfsfc,zscale,ustar,wm,ws,hbl)
!---------------------------------------------------------------------------------------------------
      implicit none
      real(8),intent(in)  :: hbl
      real(8)             ::  Bfsfc,zscale,ustar, wm,ws
      real(8)             ::  ustar3,zetahat,zeta_m,zeta_s
      real(8)             ::  a_m,c_m,a_s,c_s,r2,r3,r4
      real(8), parameter :: epssfc =  0.1
      real(8), parameter :: vonKar = 0.4

      parameter (                              &
          zeta_m=-0.2,  &  ! Maximum stability parameters "zeta"
          zeta_s=-1.0,  &  ! value of the 1/3 power law regime of
                          ! flux profile for momentum and tracers
          a_m=1.257,    &
          a_s=-28.86,   &   ! Coefficients of flux profile
          c_m=8.360,    &   ! for momentum and tracers in their
          c_s=98.96,    &   ! 1/3 power law regime;
          r2=0.5, r3=1./3., r4=0.25)
!
! Computes turbulent velocity scale for momentum and tracer
! using a 2D-lookup table as a function of "ustar" and "zetahat".
!
! Input:  Bfsfc
!         zscale   boundary layer depth [m].
!         ustar
!
! Output: wm,ws  turbulent velocity scales [m/s] at zscale
!                for momentum and tracer fields respectively.
!
! This routine was adapted from Bill Large 1995 code.
!-------
            if (Bfsfc.lt.0.) zscale=min(zscale, hbl*epssfc)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar**3
!
! Stable regime
!-------
            if (zetahat .ge. 0.) then
              wm=vonKar*ustar*ustar3/max( ustar3+5.*zetahat,1.E-20 )
              ws=wm
            else
!
! Unstable regime
!-------
              if (zetahat .gt. zeta_m*ustar3) then
                wm=vonKar*( ustar*(ustar3-16.*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat .gt. zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16.*zetahat)/ustar )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
!---------------------------------------------------------------------------------------------------
      end subroutine lmd_wscale
!===================================================================================================


!===================================================================================================
subroutine alfabeta(N,alpha,beta,t,rho0)
!---------------------------------------------------------------------------------------------------
      integer, intent(in   )         :: N
      real(8), intent(in   )         :: t(1:N,2),rho0
      real(8), intent(  out)         :: alpha,beta
! local variables
      real(8)    ZQ01, ZQ02, ZQ03, ZQ04, ZQ05, ZU00, ZU01, ZU02, ZU03,     &
                ZU04, ZV00, ZV01, ZV02, ZW00,    Tt, Ts, sqrtTs, cff
      parameter(ZQ01=+6.793952E-2, ZQ02=-9.095290E-3,  &
                ZQ03=+1.001685E-4, ZQ04=-1.120083E-6, ZQ05=+6.536332E-9,  &
                ZU00=+0.824493   , ZU01=-4.08990E-3 , ZU02=+7.64380E-5 ,  &
                ZU03=-8.24670E-7 , ZU04=+5.38750E-9 , ZV00=-5.72466E-3 ,  &
                ZV01=+1.02270E-4 , ZV02=-1.65460E-6 , ZW00=+4.8314E-4  )
!---------------------------------------------------------------------------------------------------
      Tt=t(N,1)
      Ts=t(N,2)
      sqrtTs=sqrt(Ts)
      cff=1./rho0
      alpha=-cff*(ZQ01+Tt*(2.*ZQ02+Tt*(3.*ZQ03+Tt*(4.*ZQ04+Tt*5.*ZQ05)))   &
                          +Ts*(ZU01+Tt*(2.*ZU02+Tt*(3.*ZU03+Tt*4.*ZU04))  &
                                            +sqrtTs*(ZV01+Tt*2.*ZV02))  &
                                                                  )
      beta= cff*(ZU00+Tt*(ZU01+Tt*(ZU02+Tt*(ZU03+Tt*ZU04)))+1.5*(ZV00+      &
                              Tt*(ZV01+Tt*ZV02) )*sqrtTs+2.*ZW00*Ts)
!---------------------------------------------------------------------------------------------------
end subroutine alfabeta
!===================================================================================================




!===================================================================================================
      SUBROUTINE lmd_bkpp(N,Akv,Akt,hbls,u,v,bvf,z_r,z_w,Hz,Ricr,f,r_D,Zob)
!---------------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in   )                   :: N
      real(8), intent(inout)                   :: Akv     (0:N  )
      real(8), intent(inout)                   :: Akt     (0:N,2)
      real(8), intent(in   )                   :: u       (1:N  )
      real(8), intent(in   )                   :: v       (1:N  )
      real(8), intent(in   )                   :: bvf     (0:N  )
      real(8), intent(in   )                   :: z_r     (1:N  )
      real(8), intent(in   )                   :: Hz      (1:N  )
      real(8), intent(in   )                   :: z_w     (0:N  )
      real(8), intent(in   )                   :: hbls
      real(8), intent(in   )                   :: f,r_D,Zob,Ricr
! local variables
      integer                                  :: k,kbl
      real(8)     :: Cr(0:N),zsbl,zscale,hbbl,zbl
      real(8)     :: ws,wm,sigma,Kern,cff,Kv0,Kt0,Ks0
      real(8)            :: ustar_bot
!      real(8), parameter :: Ricr   =  0.45
      real(8)            :: Ri_inv
      real(8), parameter :: C_Ek   = 258.
      real(8), parameter :: eps    = 1.D-20
      real(8), parameter :: g      =  9.81
      real(8), parameter :: vonKar =  0.4
       real(8)     :: Kv_bl,dKv_bl,Gm1,dGm1dS,cff_up,cff_dn
       real(8)     :: Kt_bl,dKt_bl,Gt1,dGt1dS,lmd_a1,lmd_a2,lmd_a3
       real(8)     :: Ks_bl,dKs_bl,Gs1,dGs1dS,Gm,Gs,Gt,sig
!---------------------------------------------------------------------------------------------------
      Ri_inv    =  1./Ricr
      ustar_bot = sqrt( r_D * sqrt( u(1)**2 + v(1)**2  ) )
      wm        = vonKar*ustar_bot
      ws        = wm
      kbl       = N
      Cr(0)     = 0.
!
!-----------------------------------------------------------------
!  Compute the Integral
!-----------------------------------------------------------------
!
      DO k=1,N-1,+1
         zscale=z_r(k)-z_w(0)
         Kern=zscale/(zscale+Zob)
         Cr(k)=Cr(k-1) + Kern*(                                                                   &
              2.0*( ( u(k+1)-u(k) )**2 + ( v(k+1)-v(k) )**2 )/(Hz(k)+Hz(k+1)) &
             -0.5*(Hz(k)+Hz(k+1))*( Ri_inv*bvf(k)+C_Ek*f*f ))
      ENDDO

      Cr(N)=2.*Cr(N-1)-Cr(N-2)  ! extrapolation toward the bottom

      DO k=1,N,+1
         if (kbl.eq.N .and. Cr(k).lt.0.) kbl=k
      ENDDO
!
!-----------------------------------------------------------------
! Linear interpolation to find hbl
!-----------------------------------------------------------------
!

      hbbl = z_w(N)-z_w(0)
      IF (kbl.lt.N) THEN
         k=kbl
         IF (k.eq.1) THEN
            hbbl=z_r(1)-z_w(0)
         ELSE
            hbbl=( z_r(k)*Cr(k-1)-z_r(k-1)*Cr(k) )/(Cr(k-1)-Cr(k)) -z_w(0)
         END IF
      END IF
!
!-----------------------------------------------------------------
! Compute boundary layer mixing coefficients.
!-----------------------------------------------------------------
!
!      DO k=1,N-1
!         IF (k.lt.kbl) THEN
!            sigma=(z_w(k)-z_w(0))/(hbbl+eps)
!            IF (sigma.lt.1.) then
!               cff=sigma*(1.-sigma)**2
!            ELSE
!               cff=0.
!            ENDIF
!            Kv0=cff*wm*hbbl
!            Kt0=cff*ws*hbbl
!            Ks0=cff*ws*hbbl
!            zsbl=z_w(N)-hbls
!            IF (z_w(k).gt.zsbl) THEN
!               Kv0=max(AKv(k  ),Kv0)
!               Kt0=max(AKt(k,1),Kt0)
!               Ks0=max(AKt(k,2),Ks0)
!            ENDIF
!            AKv(  k)=Kv0
!            AKt(k,1)=Kt0
!            AKt(k,2)=Ks0
!         ENDIF
!      ENDDO


      zbl=z_w(0)+hbbl
      k=kbl
      if (zbl.lt.z_w(k-1)) k=k-1
      cff=1./(z_w(k)-z_w(k-1))
      cff_up=cff*(zbl-z_w(k-1))
      cff_dn=cff*(z_w(k)-zbl)

      Kv_bl=cff_up*AKv(k)+cff_dn*AKv(k-1)
      dKv_bl=-cff*(AKv(k)-AKv(k-1))
      Gm1=Kv_bl/(hbbl*wm+eps)
      dGm1dS=min(0.,dKv_bl/(wm+eps))

      Kt_bl=cff_up*AKt(k,1)+cff_dn*AKt(k-1,1)
      dKt_bl=-cff*(AKt(k,1)-AKt(k-1,1))
      Gt1=Kt_bl/(hbbl*ws+eps)
      dGt1dS=min(0.,dKt_bl/(ws+eps))

      Ks_bl=cff_up*AKt(k,2)+cff_dn*AKt(k-1,2)
      dKs_bl=-cff*(AKt(k,2)-AKt(k-1,2))
      Gs1=Ks_bl/(hbbl*ws+eps)
      dGs1dS=min(0.,dKs_bl/(ws+eps))


!
!-----------------------------------------------------------------
!  Compute boundary layer mixing coefficients.
!-----------------------------------------------------------------
!
      do k=1,N-1

         if (k.lt.kbl) then
!
!  Set polynomial coefficients for shape function.
!
            sig=min((z_w(k)-z_w(0))/(hbbl+eps),1.)
            lmd_a1=sig-2.
            lmd_a2=3.-2.*sig
            lmd_a3=sig-1.
!
!  Compute nondimensional shape functions.
!
            Gm=lmd_a1+lmd_a2*Gm1+lmd_a3*dGm1dS
            Gt=lmd_a1+lmd_a2*Gt1+lmd_a3*dGt1dS
            Gs=lmd_a1+lmd_a2*Gs1+lmd_a3*dGs1dS
!
!  Compute boundary layer mixing coefficients, combine them
!  with interior mixing coefficients.
!
            Kv0=hbbl*wm*sig*(1.+sig*Gm)
            Kt0=hbbl*ws*sig*(1.+sig*Gt)
            Ks0=hbbl*ws*sig*(1.+sig*Gs)

            zsbl=z_w(N)-hbls

            if (z_w(k).gt.zsbl) then
               Kv0=max(Kv0,AKv(k  ))
               Kt0=max(Kt0,AKt(k,1))
               Ks0=max(Ks0,AKt(k,2))
            endif

            AKv(k  )=Kv0
            AKt(k,1)=Kt0
            AKt(k,2)=Ks0

         endif

      enddo ! <-- k

!---------------------------------------------------------------------------------------------------
end subroutine lmd_bkpp
!===================================================================================================




























!===================================================================================================
subroutine lmd_kpp_LMD94(N,Akv,Akt,hbl,u,v,t,bvf,rho1,z_r,z_w,Hz,Ricr,f,sustr,svstr,srflx,stflx,rho0,swr_frac,ghat)        !<-- formulation implemented in Roms-Agrif
!---------------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in   )                   :: N
      real(8), intent(  out)                   :: Akv     (0:N  )
      real(8), intent(  out)                   :: Akt     (0:N,2)
      real(8), intent(in   )                   :: u       (1:N  )
      real(8), intent(in   )                   :: v       (1:N  )
      real(8), intent(in   )                   :: t       (1:N,2)
      real(8), intent(in   )                   :: bvf     (0:N  )
      real(8), intent(in   )                   :: z_r     (1:N  )
      real(8), intent(in   )                   :: Hz      (1:N  )
      real(8), intent(in   )                   :: rho1    (1:N  )
      real(8), intent(in   )                   :: z_w     (0:N  )
      real(8), intent(in   )                   :: swr_frac(0:N  )
      real(8), intent(inout)                   :: ghat    (0:N  )
      real(8), intent(inout)                   :: hbl
      real(8), intent(in   )                   :: Ricr,f,rho0
      real(8), intent(in   )                   :: sustr,svstr
      real(8), intent(in   )                   :: srflx,stflx(2)




       integer     :: k,kbl,ka,ku,ksave
       real(8)     :: Cg,Vtc,Rib(2),Vtsq,zscale
       real(8)     :: hekman,hmonob,Bfsfc,Bo,Bosol
       real(8)     :: ws,wm,cff_up,cff_dn,sigma,dVsq,zbl
       real(8)     :: a1,a2,a3,alpha,beta,Ritop,cff,ustar
       real(8)     :: Kv_bl,dKv_bl,Gm1,dGm1dS
       real(8)     :: Kt_bl,dKt_bl,Gt1,dGt1dS
       real(8)     :: Ks_bl,dKs_bl,Gs1,dGs1dS
       real(8), parameter :: epssfc =  0.1
       real(8), parameter :: betaT  = -0.2
       real(8), parameter :: c_s    = 98.96
       real(8), parameter :: CEkman = 0.7
       real(8), parameter :: CMonob = 1.0
       real(8), parameter :: Cstar  = 10.
       real(8), parameter :: Cv     = 1.8
       real(8), parameter :: vonKar = 0.4
       real(8), parameter :: g      = 9.81
       real(8), parameter :: eps    = 1.E-14
!
! Nondimensional constants for computing non-local flux and convective
! deepening of surface boundary layer.
!------

      Cg =Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT)/( sqrt(c_s*epssfc)*Ricr*vonKar*vonKar)
!
! compute turbulent  friction velocity "ustar" from wind stress
!------

      ustar=sqrt( sqrt(sustr**2+svstr**2) )
!
! Compute thermal expansion, "alpha" [kg/m^3/degC], and saline
! contraction, "beta" [kg/m^3/PSU], coefficients at surface
!------

      call alfabeta(N,alpha,beta,t,rho0)

!
! compute surface turbulent buoyancy forcing "Bo" [m^2/s^3] (in doing
! so remove incoming solar shortwave radiation component and save it
! separately as "Bosol")
!------
      Bo   =g*(alpha*(stflx(1)-srflx)-beta * stflx(2))
      Bosol=g* alpha* srflx

!
!  Compute bulk Richardson number "Rib" and then find depth of the
!  oceanic planetary boundary layer "hbl", such that Rib(hbl)=Ric.
!------

      ka      = 1
      ku      = 2
      hbl     = z_w(N)-z_r(1)
      kbl     = 1
      Rib(ka) = 0.
      cff     = g/rho0

      do k=N-1,1,-1

        Bfsfc=Bo+Bosol*(1.-swr_frac(k))
        zscale = min(z_w(N)-z_r(k),1.)
        Call lmd_wscale(Bfsfc,zscale,ustar,wm,ws,hbl)

        Ritop=-cff*(rho1(N)-rho1(k))*(z_r(N)-z_r(k))

        dVsq= (u(N)-u(k))**2 +(v(N)-v(k))**2

        Vtsq= Vtc*(z_r(N)-z_r(k))*ws                &
              *sqrt(max(0.,0.5*(bvf(k)+bvf(k-1))))

        Rib(ku)=Ritop/(dVsq+Vtsq+eps)

        if (kbl.eq.1 .and. Rib(ku).gt.Ricr) then
           zbl=z_r(k+1)-(z_r(k+1)-z_r(k))*        &
                    (Ricr-Rib(ka))/(Rib(ku)-Rib(ka))
           hbl=z_w(N)-zbl
           kbl=k
        endif

        ksave=ka
        ka=ku
        ku=ksave

      enddo

!
!  Find stability and buoyancy forcing "Bfsfc" at "hbl".
!----------

      Bfsfc=Bo+Bosol*(1.-swr_frac(kbl-1))

!
!  Correct "hbl" with physically limiting cases (Ekman depth
!  and Monin-Obukhov depth).
!----------
      if (Bfsfc.gt.0.) then
         hekman=cekman*ustar/max(abs(f),eps)
         hbl   =min(hbl,hekman)

         hmonob=cmonob*ustar*ustar*ustar/(vonKar*Bfsfc)
         hbl   =min(hbl,hmonob)
      endif
!
!  Compute tubulent velocity scales (wm,ws) at "hbl".
!-----------
      sigma=hbl*epssfc
      call lmd_wscale (Bfsfc,sigma,ustar,wm,ws,hbl)

!
!  Compute nondimensional shape function Gx(sigma) at "hbl"
!  (sigma=1) in terms of interior diffusivities (Gx1) and
!  its vertical derivative (dGx1dS) via interpolation.
!-----------------------------------------------------------------
!
      zbl=z_w(N)-hbl
      k=kbl
      if (zbl.gt.z_w(k)) k=k+1
      cff   = 1./ (z_w(k)-z_w(k-1))
      cff_up= cff*(zbl   -z_w(k-1))
      cff_dn= cff*(z_w(k)-zbl)

      Kv_bl=cff_up*AKv(k)+cff_dn*AKv(k-1)
      dKv_bl=cff*(AKv(k)-AKv(k-1))

      Kt_bl=cff_up*AKt(k,1)+cff_dn*AKt(k-1,1)
      dKt_bl=cff*(AKt(k,1)-AKt(k-1,1))

      Ks_bl=cff_up*AKt(k,2)+cff_dn*AKt(k-1,2)
      dKs_bl=cff*(AKt(k,2)-AKt(k-1,2))

      Gm1=Kv_bl/(hbl*wm+eps)
      dGm1dS=min(0.,-dKv_bl/(wm+eps))

      Gt1=Kt_bl/(hbl*ws+eps)
      dGt1dS=min(0.,-dKt_bl/(ws+eps))

      Gs1=Ks_bl/(hbl*ws+eps)
      dGs1dS=min(0.,-dKs_bl/(ws+eps))

! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.

       ghat(0:N) = 0.

       do k=N-1,kbl,-1

         sigma=z_w(N)-z_w(k)
         Call lmd_wscale(Bfsfc,sigma,ustar,wm,ws,hbl)
! Compute vertical mixing coefficients
         sigma=(z_w(N)-z_w(k))/max(hbl,eps)
         a1=sigma-2.
         a2=3.-2.*sigma
         a3=sigma-1.

         AKv(k)=wm*hbl*(       sigma*( 1.+sigma*(    &
                                 a1+a2*Gm1+a3*dGm1dS )))
         AKt(k,1)=ws*hbl*(       sigma*( 1.+sigma*(    &
                                 a1+a2*Gt1+a3*dGt1dS )))
         AKt(k,2)=ws*hbl*(       sigma*( 1.+sigma*(    &
                                 a1+a2*Gs1+a3*dGs1dS )))
! Compute non-local term
         if (Bfsfc .lt. 0.) then
            !ghat(k)= Cg * sigma*( 1.+sigma*( a1+a2*Gt1+a3*dGt1dS ) )
            ghat(k)=Cg * sigma*(1.-sigma)**2
         else
            ghat(k)=0.
         endif

       enddo

!---------------------------------------------------------------------------------------------------
end subroutine lmd_kpp_LMD94
!===================================================================================================







!===================================================================================================
      SUBROUTINE lmd_bkpp_LMD94(N,Akv,Akt,hbl,u,v,bvf,rho0,rho1,z_r,z_w,Hz,Ric,f,r_D,Zob)
!---------------------------------------------------------------------------------------------------
      implicit none
      integer, intent(in   )                   :: N
      real(8), intent(inout)                   :: Akv     (0:N  )
      real(8), intent(inout)                   :: Akt     (0:N,2)
      real(8), intent(in   )                   :: u       (1:N  )
      real(8), intent(in   )                   :: v       (1:N  )
      real(8), intent(in   )                   :: bvf     (0:N  )
      real(8), intent(in   )                   :: z_r     (1:N  )
      real(8), intent(in   )                   :: Hz      (1:N  )
      real(8), intent(in   )                   :: rho1    (1:N  )
      real(8), intent(in   )                   :: z_w     (0:N  )
      real(8), intent(in   )                   :: hbl,rho0
      real(8), intent(in   )                   :: f,r_D,Zob,Ric
! local variables
      integer                                  :: k,kbbl,ka,ku,ksave
      real(8)            :: Rib(2),zsbl,zscale,hbbl,Ritop,zbl
      real(8)            :: ws,wm,sig,cff,Kv0,Kt0,Ks0
      real(8)            :: dGm1dS,Gm1,Gt1,dGt1dS,Gs1,dGs1dS,Kv_bl,Kt_bl,Ks_bl
      real(8)            :: dKv_bl,dKt_bl,dKs_bl,Gm,Gt,Gs,lmd_a1,lmd_a2,lmd_a3
      real(8)            :: ustar ,cff_up,cff_dn,dVsq,hekman,Vtc,Vtsq
      real(8)            :: Ri_inv
      real(8), parameter :: lmd_cs=98.96
      real(8), parameter :: lmd_Cv=1.8
      real(8), parameter :: lmd_betaT=-0.2
      real(8), parameter :: lmd_epsilon=0.1
      real(8), parameter :: lmd_cekman=0.7
      real(8), parameter :: lmd_nu0c=0.1
      real(8), parameter :: eps    = 1.D-20
      real(8), parameter :: g      =  9.81
      real(8), parameter :: vonKar =  0.4
!---------------------------------------------------------------------------------------------------
      Ri_inv    =  1./Ric
      ustar     = sqrt( r_D * sqrt( u(1)**2 + v(1)**2  ) )

      Vtc=lmd_Cv*sqrt(-lmd_betaT)/( sqrt(lmd_cs*lmd_epsilon)  &
                                              *Ric*vonKar*vonKar )
      wm        = vonKar*ustar
      ws        = wm

      ka=1
      ku=2

      hbbl=z_r(N)-z_w(0)
      kbbl=N
      Rib(ku)=0.

      do k=2,N
         cff=g/rho0
         Ritop=-cff*(rho1(k)-rho1(1))*(z_r(k)-z_r(1))
         dVsq=0.25*( (u(k)-u(1)+u(k)-u(1))**2    &
                    +(v(k)-v(1)+v(k)-v(1))**2)
         Vtsq=Vtc*(z_r(k)-z_r(1))*ws*sqrt(max(0.,0.5*(bvf(k)+bvf(k-1))))
         Rib(ka)=Ritop/(dVsq+Vtsq+eps)

         if (kbbl.eq.N .and. Rib(ka).gt.Ric) then
            zbl=z_r(k)-(z_r(k)-z_r(k-1))*(Ric-Rib(ka))/(Rib(ku)-Rib(ka))
            hbbl=zbl-z_w(0)
            kbbl=k
         endif
         ksave=ka
         ka=ku
         ku=ksave
      enddo

      hekman= lmd_cekman*ustar/max(abs(f),1.e-6)
      hbbl  = min(hbbl,hekman)

      zbl=z_w(0)+hbbl
      k=kbbl
      if (zbl.lt.z_w(k-1)) k=k-1
      cff=1./(z_w(k)-z_w(k-1))
      cff_up=cff*(zbl-z_w(k-1))
      cff_dn=cff*(z_w(k)-zbl)

      Kv_bl=cff_up*AKv(k)+cff_dn*AKv(k-1)
      dKv_bl=-cff*(AKv(k)-AKv(k-1))
      Gm1=Kv_bl/(hbbl*wm+eps)
      dGm1dS=min(0.,dKv_bl/(wm+eps))

      Kt_bl=cff_up*AKt(k,1)+cff_dn*AKt(k-1,1)
      dKt_bl=-cff*(AKt(k,1)-AKt(k-1,1))
      Gt1=Kt_bl/(hbbl*ws+eps)
      dGt1dS=min(0.,dKt_bl/(ws+eps))

      Ks_bl=cff_up*AKt(k,2)+cff_dn*AKt(k-1,2)
      dKs_bl=-cff*(AKt(k,2)-AKt(k-1,2))
      Gs1=Ks_bl/(hbbl*ws+eps)
      dGs1dS=min(0.,dKs_bl/(ws+eps))


!
!-----------------------------------------------------------------
!  Compute boundary layer mixing coefficients.
!-----------------------------------------------------------------
!
      do k=1,N-1

         if (k.lt.kbbl) then
!
!  Set polynomial coefficients for shape function.
!
            sig=min((z_w(k)-z_w(0))/(hbbl+eps),1.)
            lmd_a1=sig-2.
            lmd_a2=3.-2.*sig
            lmd_a3=sig-1.
!
!  Compute nondimensional shape functions.
!
            Gm=lmd_a1+lmd_a2*Gm1+lmd_a3*dGm1dS
            Gt=lmd_a1+lmd_a2*Gt1+lmd_a3*dGt1dS
            Gs=lmd_a1+lmd_a2*Gs1+lmd_a3*dGs1dS
!
!  Compute boundary layer mixing coefficients, combine them
!  with interior mixing coefficients.
!
            Kv0=hbbl*wm*sig*(1.+sig*Gm)
            Kt0=hbbl*ws*sig*(1.+sig*Gt)
            Ks0=hbbl*ws*sig*(1.+sig*Gs)

            zsbl=z_w(N)-hbl

            if (z_w(k).gt.zsbl) then
               Kv0=max(Kv0,AKv(k  ))
               Kt0=max(Kt0,AKt(k,1))
               Ks0=max(Ks0,AKt(k,2))
            endif

            AKv(k  )=Kv0
            AKt(k,1)=Kt0
            AKt(k,2)=Ks0

         endif

      enddo ! <-- k

!---------------------------------------------------------------------------------------------------
end subroutine lmd_bkpp_LMD94
!===================================================================================================



end module scm_kpp
