"""
This module is to convert the source codes ready for umat/vumat build.
Conversely, it can be used to convert the source codes ready for stand-alone build.
"""

def sub_vaba(filename,
             init="include 'aba_param.inc'",
             fin="include 'vaba_param.inc'"):
    """
    Given the <filename>, replace the <init> phrase to
    <fin> phrase in every occurence.

    Arguments
    ---------
    filename
    init
    fin
    """
    import re
    p=re.compile(r"%s"%init,re.IGNORECASE)
    with open(filename,'r') as fo:
        cnt=fo.read()
    m=p.findall(cnt)
    s=p.sub(r"%s"%fin,cnt)
    with open(filename,'w') as fo:
        fo.write(s)

if __name__=='__main__':
    """
    python prebuild_vumat.py -i <0 or 1, -1, -2>

    Argument
    --------
    -i   : options
           i=0: aba         -> vaba
           i=1: vaba        -> aba
           i=2: SA          -> implicit FE
           i=3: SA          -> explicit FE
           i=4: implicit FE -> SA
           i=5: explicit FE -> SA
    """
    from glob import glob
    import os
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', type=int,
        help='option (i=0: aba -> vaba; \
        i=1: vaba -> aba; i=2: SA -> implicit FE; \
        i=3: SA -> explicit FE ; i=4: implicit FE -> SA; \
        i=5: explicit FE -> SA')
    args=parser.parse_args()

    #if args.i==0 or args.i==1 or args.i==2:
    ## source code files subject to be converted.
    fns=glob('*.f')+glob('src/*.f')+glob('../../../src/vpsc_src/*.f')+glob('../../../src/vpsc_src/*.sub')

    # if args.i==2:
    #      fns=fns+glob('*.f')+glob('yld/*.f')+glob('yld/hah/*.f')+glob('lib/*.f')

    for f in fns:
        if not(os.path.isfile(f)):
            print('could not find the file %s'%f)
        if args.i==0:
            sub_vaba(filename=f,init="include 'ABA_PARAM.INC'", fin="include 'VABA_PARAM.INC'")
        elif args.i==1:
            sub_vaba(filename=f,init="include 'VABA_PARAM.INC'",fin="include 'ABA_PARAM.INC'" )
        elif args.i==2:
            sub_vaba(filename=f,init="implicit none",fin="include 'ABA_PARAM.INC'" )
        elif args.i==3:
            sub_vaba(filename=f,init="implicit none",fin="include 'VABA_PARAM.INC'" )
        elif args.i==4:
            sub_vaba(filename=f,fin="implicit none",init="include 'ABA_PARAM.INC'" )
        elif args.i==5:
            sub_vaba(filename=f,fin="implicit none",init="include 'VABA_PARAM.INC'" )
