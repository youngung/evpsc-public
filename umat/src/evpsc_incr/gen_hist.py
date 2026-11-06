## python script to generate velocity gradient history for stand-alone EVPSC calculation.
## For details, take a look at its companion Python notebook: <gen_hist.ipynb>

## Python 3 compatible.

import numpy as np
import matplotlib.pyplot as plt

def readnew(fn='tmp_0001_0001.rst',ntens=6):
    """
    Arguments
    ---------
    fn     ='tmp_0001_0001.rst'
    ntens  =6

    Returns
    -------
    udot
    """
    import sys
    print('** starts readnew')
    print('fn:',fn)
    dat=np.loadtxt(fn,skiprows=1).T
    time=dat[2,:] # time
    dat=dat[:,1:] ## remove the first row
    dt=np.diff(time) # time increments
    print(f'initial {dt.shape},{dat.shape}')

    ## remove data with zero time increment.
    filt=dt>1e-24
    dat=dat[:,filt]

    ## number of increments
    dt=dt[filt]
    nincrs=dat.shape[-1]
    print(f'later {dt.shape},{dat.shape}')

    ## strip off data columns that are precendent to UVARM
    i0=6
    i1=i0+ntens
    print('i0, i1:',i0,i1)
    LE=dat[i0:i1]
    i0=i1
    i1=i0+ntens
    print('i0, i1:',i0,i1)
    sig=dat[i0:i1]
    i0=i1
    i1=i0+ntens
    print('i0, i1:',i0,i1)
    edot=dat[i0:i1]
    uvarm=dat[i1:]

    udot=uvarm[:9,:] ## velocity gradient

    records=np.zeros((10,nincrs))
    records[0:9,:]=udot[:]
    records[9,:]=dt[:]
    return records.T


def gennew(records,fnhist='/tmp/lij.hist'):
    """
    Arguments
    ---------
    records
    fnhist
    """
    import sys
    #with open(fnhist,'w') as fo:
    if True:
        print('records.shape:',records.shape)
        nincr,ncol=records.shape
        cnt='%5s'%'step'
        items=['l11','l12','l13','l21','l22','l23','l31','l32','l33','tincr']
        for item in items:
            cnt=f'{cnt}%20s'%item
        cnt=f"{cnt}\n"

        for i in range(nincr):
            cnt=f"{cnt}%5i"%(i+1)
            for j in range(9):
                cnt=f"{cnt}%20.12e"%records[i,j]
            cnt=f"{cnt}%20.12e\n"%records[i,9]

        cnt=f'{nincr} 7 0 298.\n{cnt}'

        print('---')
        print(cnt)
        print('---')

        with open(fnhist,'w') as fo:
            fo.write(cnt)

    print(f"{fnhist} has been created")

def main(fnvel,fnhist,ntens):
    records = readnew(fnvel,ntens)
    gennew(records,fnhist)

if __name__=='__main__':
    import argparse, sys
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',type=str,help='Result file extracted from Abaqus ODB file generated from FE EVPSC simulations (UMAT)',default='tmp_0001_0001.rst')
    parser.add_argument('-o',type=str,help='History file to be generated',default='/tmp/lij.txt')
    parser.add_argument('-ntens',type=int,help='ntens')

    args = parser.parse_args()

    if type(args.ntens)==type(None):
        msg='** Error: argument <ntens> is missing in <gen_hist.py>'
        sys.stderr.write(f'{msg}\n\n')
        raise IOError

    main(args.i,args.o,args.ntens)
