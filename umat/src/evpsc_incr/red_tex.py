def mmm_symop():
    import numpy as np
    m0 = [[ 1, 0, 0], [0, 1, 0], [0, 0, 1]]
    m1 = [[ 1, 0, 0], [0,-1, 0], [0, 0,-1]]
    m2 = [[-1, 0, 0], [0, 1, 0], [0, 0,-1]]
    m3 = [[-1, 0, 0], [0,-1, 0], [0, 0, 1]]
    h = np.array([m0,m1,m2,m3])
    return h



def euler(ph=None, th=None, tm=None, a=None, echo=True):
    """
    note:
          This is a pythonized fortran subroutine embedded in the VPSC7.sub
          done by youngung 2010-12-23

          Matrix A(i,j) is returned as a numpy array
          all Euler angles are in degree not in radian

          # A[i,j] transforms from SA to CA

          Thanks to python's non-hardwirable-arugment feature,
          if a matrix is not given, it automatically calculates
          the matrix with given angles.
          Vice versa, if matrix is given, given angle aurgments
          are ignored and new euler angles are returned.

    Nomenclature of Euler angle follows Bunge's convention
          ph = phi1,
          th = phi
          tm = phi2
    """
    import numpy as np
    import math
    if type(a).__name__=='NoneType':  a=np.resize(np.array(()),(3,3));iopt=2
    else:
        if type(a).__name__=='ndarray':
            iopt = 1
            pass
        else:
            print('Error: Unexpected Matrix a type')
            print('It should be numpy.ndarry!')
            raise IOError

    if iopt==1 :
        th = math.acos(a[2,2])  #Radian
        if abs(a[2,2] > 0.99999):
            tm = 0.
            ph = math.atan2(a[0,1],a[0,0]) #Radian
        else:
            sth = math.sin(th)
            tm = math.atan2(a[0,2]/sth,a[1,2]/sth)
            ph = math.atan2(a[2,0]/sth,-a[2,1]/sth)
        th = th * 180./math.pi
        ph = ph * 180./math.pi
        tm = tm * 180./math.pi
        return [ph,th,tm] #phi1, phi, phi2

    elif (iopt == 2):
        angles = [ph,th,tm]
        if any(angles[i] == None for i in range(len(angles))):
            print('Angles must be give if iopt==2')
            raise IOError

        """ Convert the angle into Radian"""
        ph = ph * math.pi / 180.
        th = th * math.pi / 180.
        tm = tm * math.pi / 180.

        sph = math.sin(ph)
        cph = math.cos(ph)
        sth = math.sin(th)
        cth = math.cos(th)
        stm = math.sin(tm)
        ctm = math.cos(tm)

        a[0,0] =  ctm * cph - sph * stm * cth
        a[1,0] = -stm * cph - sph * ctm * cth
        a[2,0] =  sph * sth
        a[0,1] =  ctm * sph + cph * stm * cth
        a[1,1] = -sph * stm + cph * ctm * cth
        a[2,1] = -sth * cph
        a[0,2] =  sth * stm
        a[1,2] =  ctm * sth
        a[2,2] =  cth

        if echo==True:
            print('Matrix a is ')
            print(a)
        return a

    else: print('Error: Improper iopt input'); raise IOError


if __name__=='__main__':
    import numpy as np
    import tempfile
    import argparse
    fn=tempfile.mktemp(prefix="youngung_tmp_")
    # print('fn:',fn)
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', type=str, help='input file')
    parser.add_argument(
        '-o', type=str, help='output file',default=fn)
    parser.add_argument(
        '-n', type=int, help='Number of grains',default=1)
    parser.add_argument(
        '-s', action='store_true', help='Apply orthotropic symmetry')

    args = parser.parse_args()

    try:
        grs_old=np.loadtxt(args.i,skiprows=4)
    except:
        with open(args.i,'r') as fo:
            lines=fo.read().split('\n')
            ngr=int(lines[3].split()[1])
            lines=lines[4:4+ngr]
            grs_old=np.zeros((ngr,4))
            for i, line in enumerate(lines):
                grs_old[i,:]=np.array(line.split(),dtype='float')

    number_of_rows = grs_old.shape[0]

    if args.s:
        remainder=np.mod(args.n,4)
        # print('remainder:',remainder)
        # print('args.n:',args.n)
        # print('number_of_rows:',number_of_rows)
        if remainder!=0:
            msg=''
            msg=f'{msg}** When applying the orthotropic symmetry,\n'
            msg=f'{msg}** the total number of requested orientations (n)\n'
            msg=f'{msg}** should be a multiple of 4'
            raise IOError(msg)
        args.n=int(args.n/4.)

    random_indices = np.random.choice(number_of_rows, size=args.n, replace=False)


    if args.s:
        mmm=mmm_symop()
        grs_ref = grs_old[random_indices,:]
        grs_new = np.zeros((args.n*4,4))
        for i, ori in enumerate(grs_ref):
            cRs=euler(*ori[:3],echo=False)
            for j, mop in enumerate(mmm):
                R0=np.dot(cRs, mop)
                phi1, phi, phi2 = euler(a=R0, echo=False)
                grs_new[i*4+j,:]=phi1,phi,phi2,ori[3]
    else:
        grs_new = grs_old[random_indices, :]


    wgt=grs_new[:,3]
    wgt=wgt/wgt.sum()
    grs_new[:,3]=wgt[:]

    if args.s:
        args.n=int(args.n*4)

    with open(args.o,'w') as fo:
        fo.write("1\n2\n3\nB %i\n"%args.n)
        for i in range(args.n):
            fo.write("%8.3f %8.3f %8.3f %16.8e \n"%(
                grs_new[i,0],grs_new[i,1],grs_new[i,2],grs_new[i,3]))
