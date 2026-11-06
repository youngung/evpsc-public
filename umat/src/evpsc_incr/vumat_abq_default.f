       subroutine vumat(
C Read only -
     $     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $     stepTime, totalTime, dt, cmname, coordMp, charLength,
     $     props, density, strainInc, relSpinInc,
     $     tempOld, stretchOld, defgradOld, fieldOld,
     $     stressOld, stateOld, enerInternOld, enerInelasOld,
     $     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     $     stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
C J2 Mises Plasticity with kinematic hardening for plane
C strain case.
C Elastic predictor, radial corrector algorithm.
C
C The state variables are stored as:
C      STATE(*,1) = back stress component 11
C      STATE(*,2) = back stress component 22
C      STATE(*,3) = back stress component 33
C      STATE(*,4) = back stress component 12
C      STATE(*,5) = equivalent plastic strain
C
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     $     coordMp(nblock,*),
     $     charLength(*), strainInc(nblock,ndir+nshr),
     $     relSpinInc(*), tempOld(*),
     $     stretchOld(*), defgradOld(*),
     $     fieldOld(*), stressOld(nblock,ndir+nshr),
     $     stateOld(nblock,nstatev), enerInternOld(nblock),
     $     enerInelasOld(nblock), tempNew(*),
     $     stretchNew(*), defgradNew(*), fieldNew(*),
     $     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     $     enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
      parameter( zero = 0., one = 1., two = 2., three = 3.,
     $     third = one/three, half = .5, twoThirds = two/three,
     $     threeHalfs = 1.5 )
C
      e     = props(1)
      xnu   = props(2)
      yield = props(3)
      hard  = props(4)
C
      twomu  = e / ( one + xnu )
      thremu = threeHalfs * twomu
      sixmu  = three * twomu
      alamda = twomu * ( e - twomu ) / ( sixmu - two * e )
      term   = one / ( twomu * ( one + hard/thremu ) )
      con1   = sqrt( twoThirds )
C
      do 100 i = 1,nblock
C
C Trial stress
        trace  = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
        sig1 = stressOld(i,1) + alamda*trace + twomu*strainInc(i,1)
        sig2 = stressOld(i,2) + alamda*trace + twomu*strainInc(i,2)
        sig3 = stressOld(i,3) + alamda*trace + twomu*strainInc(i,3)
        sig4 = stressOld(i,4)                + twomu*strainInc(i,4)
C
C Trial stress measured from the back stress
        s1 = sig1 - stateOld(i,1)
        s2 = sig2 - stateOld(i,2)
        s3 = sig3 - stateOld(i,3)
        s4 = sig4 - stateOld(i,4)
C
C Deviatoric part of trial stress measured from the back stress
        smean = third * ( s1 + s2 + s3 )
        ds1 = s1 - smean
        ds2 = s2 - smean
        ds3 = s3 - smean
C
C Magnitude of the deviatoric trial stress difference
        dsmag = sqrt( ds1**2 + ds2**2 + ds3**2 + 2.*s4**2 )
C
C Check for yield by determining the factor for plasticity,
C zero for elastic, one for yield
        radius = con1 * yield
        facyld = zero
        if(  dsmag - radius .ge. zero ) facyld = one
C
C Add a protective addition factor to prevent a divide by zero
C when dsmag is zero. If dsmag is zero, we will not have exceeded
C the yield stress and facyld will be zero.
        dsmag  = dsmag + ( one - facyld )
C
C Calculated increment in gamma (this explicitly includes the
C time step)
        diff   = dsmag - radius
        dgamma = facyld * term * diff
C
C Update equivalent plastic strain
        deqps  = con1 * dgamma
        stateNew(i,5) = stateOld(i,5) + deqps
C
C Divide dgamma by dsmag so that the deviatoric stresses are
C explicitly converted to tensors of unit magnitude in the
C following calculations
        dgamma = dgamma / dsmag
C
C Update back stress
        factor  = hard * dgamma * twoThirds
        stateNew(i,1) = stateOld(i,1) + factor * ds1
        stateNew(i,2) = stateOld(i,2) + factor * ds2
        stateNew(i,3) = stateOld(i,3) + factor * ds3
        stateNew(i,4) = stateOld(i,4) + factor *  s4
C
C Update the stress
        factor   = twomu * dgamma
        stressNew(i,1) = sig1 - factor * ds1
        stressNew(i,2) = sig2 - factor * ds2
        stressNew(i,3) = sig3 - factor * ds3
        stressNew(i,4) = sig4 - factor *  s4
C
C Update the specific internal energy -
        stressPower = half * (
     $    ( stressOld(i,1)+stressNew(i,1) )*strainInc(i,1)
     $    +     ( stressOld(i,2)+stressNew(i,2) )*strainInc(i,2)
     $    +     ( stressOld(i,3)+stressNew(i,3) )*strainInc(i,3)
     $    + two*( stressOld(i,4)+stressNew(i,4) )*strainInc(i,4) )
C
        enerInternNew(i) = enerInternOld(i)
     $    + stressPower / density(i)
C
C Update the dissipated inelastic specific energy -
        smean = third * ( stressNew(i,1) + stressNew(i,2)
     $                + stressNew(i,3) )
        equivStress = sqrt( threeHalfs *
     $    ( (stressNew(i,1)-smean)**2
     $    + (stressNew(i,2)-smean)**2
     $    + (stressNew(i,3)-smean)**2
     $    + two * stressNew(i,4)**2 ) )
        plasticWorkInc = equivStress * deqps
        enerInelasNew(i) = enerInelasOld(i)
     $       + plasticWorkInc / density(i)
  100 continue
C
      return
      end
