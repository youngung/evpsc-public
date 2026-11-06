"""
Generate dimension file.
Usage.
python -u pyscripts/gen_dim.py -o src/parameters.dim --nphmx 1 --ngrmx 10 --ntwmmx 1 --nsysmx 13 --nneignmx 1 --ndifmx 1
"""
def read_fnsx(fnsx,v):
    """
    Arguments
    ---------

    Returns
    -------
    returns nmodmx, ntwmmx, nsysmx, ntwsmx
    """
    import sys
    import subprocess


    if sys.platform=='darwin':
        fnexe='sx_parser'
    elif sys.platform=='win32':
        fnexe='sx_parser.exe'

    # if sx_parser is not found, make one.
    if os.path.isfile(fnexe):  pass
    else:
        print('** making parsers')
        p=subprocess.Popen(
            'make parsers'.split(),
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdo,stde=p.communicate()
        p.wait()
        if p.poll()!=0:
            print('** Error: Building parser failed')
            print(stdo.decode())
            print(stde.decode())
            sys.exit(-1)

    if sys.platform=='darwin':
        cmd=f'./sx_parser {fnsx}'
    elif sys.platform=='win32':
        cmd=f'sx_parser {fnsx}'

    if v:
        print(f'** command used to parse the given fnsx: {fnsx}')
        print(f'cmd: {cmd}')
        print(f'sys.platform: {sys.platform}')
        print(f'os.getcwd(): {os.getcwd()}')

    p=subprocess.Popen(
        cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    stdo,stde=p.communicate()
    if v or True:
        print(f'** stdo and stde from cmd {cmd}')
        print(stdo.decode())
        print(stde.decode())

    lines=stdo.decode().split('\n')
    lines=lines[::-1]
    lines=lines[:7] ## 6, but the last one has empty..
    lines=lines[::-1]
    if v: print('lines:',lines)
    rst=dict(name='dum')
    if v: print('reading keys...')

    nmodmx=int(lines[0].split(':')[1])
    ntwmmx=int(lines[1].split(':')[1])
    nsysmx=int(lines[2].split(':')[1])
    ntwsmx=int(lines[3].split(':')[1])
    ihardlaw=int(lines[4].split(':')[1])
    isx_ver=int(lines[5].split(':')[1])

    return nmodmx, ntwmmx, nsysmx, ntwsmx, ihardlaw

if __name__=='__main__':
    import argparse
    import os
    import sys
    import shutil
    import tempfile
    ## Arguments parsing
    parser = argparse.ArgumentParser()
    tmpo=tempfile.mktemp(prefix="youngung_tmp_")
    parser.add_argument('-o',type=str,default=tmpo,help='path-to-output-file')
    #
    parser.add_argument('--fnsx',nargs='*',default=None,help='Name(s) of single crystal file(s)')
    #
    parser.add_argument('--nphmx',type=int,default=1,help='Max no. of phases')
    parser.add_argument('--ngrmx',type=int,default=10000,help='Max no. of grains')

    # parameters that can be read from single crystal file.
    parser.add_argument('--nmodmx',type=int,default=6,help='Max no. of slip and twin modes')
    parser.add_argument('--ntwmmx',type=int,default=2,help='Max no. of twin modes')
    parser.add_argument('--nsysmx',type=int,default=48,help='Max no. of slip+twin systems')
    parser.add_argument('--ntwsmx',type=int,default=24,help='Max no. of twin systems')

    # non-essential arguments
    parser.add_argument('--nneimx',type=int,default=1,help='Max no. of neighbors')
    parser.add_argument('--ndifmx',type=int,default=12,help='Max no. of diffracting planes')

    # switches/options for numerical conditions
    parser.add_argument('--iregime',type=int,default=2,help='1: VPSC: 2: delta-EVPSC; 3: sigma-EVPSC')

    ## etcs.
    parser.add_argument('--v',action='store_true',help='print out the resulting dimension file')
    args        = parser.parse_args()

    if type(args.fnsx)!=type(None):
        args.nphmx =0
        args.nmodmx=0 ## this is max per phase.
        args.ntwmmx=0 ## this is max per phase
        args.nsysmx=0 ## this is max per phase
        args.ntwsmx=0 ## this is max per phase
        for i, fnsx in enumerate(args.fnsx):
            args.nphmx =args.nphmx+1

            if not os.path.isfile(fnsx):
                print(f'**Error in <gen_dim.py>: Could not find {fnsx}')
                sys.exit(-1)

            nmodmx, ntwmmx, nsysmx, ntwsmx, ihardlaw = read_fnsx(fnsx,args.v)
            # args.nmodmx=args.nmodmx+nmodmx
            # args.ntwmmx=args.ntwmmx+ntwmmx
            # args.nsysmx=args.nsysmx+nsysmx
            # args.ntwsmx=args.ntwsmx+ntwsmx
            if nmodmx>args.nmodmx:
                args.nmodmx=nmodmx
            if ntwmmx>args.ntwmmx:
                args.ntwmmx=ntwmmx
            if nsysmx>args.nsysmx:
                args.nsysmx=nsysmx
            if ntwsmx>args.ntwsmx:
                args.ntwsmx=ntwsmx

        if args.ntwmmx==0:
            args.ntwmmx=1
        if args.ntwsmx==0:
            args.ntwsmx=1


    ##------------------------------------------------------------------
    ## Write <parameters.dim> file
    cnt=f"""c-----------------------------------------------------------------------
c     Dimension of arrays. Reduntantly large size of arrays may impact
c     the computational speed of EVPSC-UMAT.
      integer NPHMX,NGRMX,NMODMX,NTWMMX,NSYSMX,NTWSMX,NNEIMX,NDIFMX
      integer iregime
      parameter(NPHMX={args.nphmx})
      parameter(NGRMX={args.ngrmx})
      PARAMETER(NMODMX={args.nmodmx})
      PARAMETER(NTWMMX={args.ntwmmx})
      PARAMETER(NSYSMX={args.nsysmx})
      PARAMETER(NTWSMX={args.ntwsmx})
      PARAMETER(NNEIMX={args.nneimx})
      PARAMETER(NDIFMX={args.ndifmx})

c     -- 1: VPSC, 2: delta-EVPSC 3: sigma-EVPSC --
      parameter(iregime={args.iregime})
c-----------------------------------------------------------------------"""

    with open(args.o,'w') as fo:
        if args.v: print(f'writing <parameters.dim> to {args.o}')
        fo.write(cnt)
    if args.v:
        with open(args.o,'r') as fo:
            print('Resulting dimension file')
            print(fo.read())
