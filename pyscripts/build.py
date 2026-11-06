def chg_kw_in_dim(dim,kw,val):
    """
    Given the dimension file, replace the keyworded parameter (kw) to the given value (val)

    Arguments
    ---------
    dim
    kw
    val
    """
    with open(dim,'r') as fo:
        oldstring=fo.read()
        oldlines=oldstring.split('\n')

    kount=0
    for i in range(len(oldlines)):
        if 'parameter(%s='%kw in oldlines[i].lower().replace(' ',''):
            print('replacing')
            oldlines[i]=f'      parameter(%s={val})'%kw.upper()
            print('oldline:',oldlines[i])
            kount+=1

    if kount!=1:
        print('** Error: kount not is equal to 1:',kount)
        raise IOError('** Expected only one occurence for %s:%i'%(kw,kount))

    new=''
    for l in oldlines:
        new=f'{new}{l}\n'

    with open(dim,'w') as fo:
        fo.write(new)


def chg_nph_in_dim(dim, nph):
    with open(dim,'r') as fo:
        oldstring=fo.read()
        oldlines=oldstring.split('\n')

    kount=0
    for i in range(len(oldlines)):
        if 'parameter(nphmx=' in oldlines[i].lower():
            print('replacing')
            oldlines[i]=f'      parameter(nphmx={nph})'.upper()
            print('oldline:',oldlines[i])
            kount+=1

    if kount!=1:
        print('** Error: kount is greater than 1:',kount)
        raise IOError('** expected only one occurence!')

    new=''
    for l in oldlines:
        new=f'{new}{l}\n'

    with open(dim,'w') as fo:
        fo.write(new)


def chg_ngr_in_dim(dim, ngr):
    """
    change the number of grains defined in the given dimension file

    arguments
    ---------
    dim
    ngr
    """
    with open(dim,'r') as fo:
        oldstring=fo.read()
        oldlines=oldstring.split('\n')

    kount=0
    for i in range(len(oldlines)):
        if 'parameter(ngrmx=' in oldlines[i].lower():
            print('replacing')
            oldlines[i]=f'      parameter(ngrmx={ngr})'.upper()
            print('oldline:',oldlines[i])
            kount+=1

    if kount!=1:
        print('** Error: kount is greater than 1:',kount)
        raise IOError('** expected only one occurence!')

    new=''
    for l in oldlines:
        new=f'{new}{l}\n'

    with open(dim,'w') as fo:
        fo.write(new)

def chg_iodf(iodf,filename):
    import re
    print('changing iodf')
    print('iodf:',iodf)
    with open(filename,'r') as fo:
        cnt=fo.read()
    init0="iodf=0"
    init1="iodf=1"
    init2="iodf=2"
    if iodf==1:
        fin="iodf=1"
    elif iodf==2:
        fin="iodf=2"
    elif iodf==0:
        fin="iodf=0"
    else:
        raise IOError(f"** unexpected iodf : {iodf}")
    p=re.compile(r'%s|%s|%s'%(init0,init1,init2),re.IGNORECASE)
    m=p.findall(cnt); s=p.sub(r"%s"%fin,cnt)
    with open(filename,'w') as fo:
        fo.write(s)

if __name__=='__main__':
    print("""---------------------------------------------------""")
    print("""         Build stand-alone dEVPSC executable       """)
    print("""---------------------------------------------------\n\n""")

    import argparse
    import sys
    import tempfile
    import subprocess
    import shutil
    import os
    import glob
    import re

    parser = argparse.ArgumentParser()
    parser.add_argument('--iclone',action='store_true',
        help="""
        Clone and build; Warning: a cloned repo does not include the changes that are not committed.
        One should not use '--clone' to use either unstaged or uncommitted changes in the repo.

        Clone literally clones the current version and use.
        Without this option, one uses 'make zip' to do so.
        """)
    parser.add_argument('--dst',default=os.getcwd(),help='The path in which clone & build is performed; default is os.getcwd()')
    ### options associated with polycrystal characteristics.
    parser.add_argument('--dim',type=str,help='dimension file',default='src/vpsc_src/parameters.dim')
    parser.add_argument('--nphs',type=int,help='No of phases',default=1)
    parser.add_argument('--ngrs',nargs='+',default=[],help='No of grains for each file')
    ### flags.
    parser.add_argument('--ielast',action='store_true',help='make CRSS very high to suppress plasticity')
    parser.add_argument('--iodf',type=int,default=0)
    parser.add_argument('--irgvb',action='store_true',help='use RGVB')
    parser.add_argument('--idebug',action='store_true',help='Use more restrictive flags for compilation')
    parser.add_argument('--ioptfs',type=int,default=0,
                        help="0: Using the linear form of Jeong et al. (IJP, 2021) and Jeon et al. (JMRT, 2024)     1: Using the new linear form of Patra et al.")

    args=parser.parse_args()
    codec='cp949'

    if args.iclone and args.dst==os.getcwd():
        print('** Error: You cannot clone to the current working directory. Specify the path.')
        sys.exit(-1)

    if args.iclone:
        print('\n--------------------------------------------------------')
        print('** Warning!!! You are cloning the repository!!!')
        print('    Cloned repository does not include uncommitted changes')
        print('    that may exist in the original repository.')
        print('-----------------------------------------------------------\n')
        print(f'* dst: {args.dst}')
        if os.path.isdir(args.dst):
            print('*** Warning!!! The given <dst> already exists')
            print('*** I will be removing the existing directory')
            shutil.rmtree(args.dst)

        cmd=f'git clone {os.getcwd()} {args.dst}'
        p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdo,stde=p.communicate()
        # print(stdo.decode())
        # print(stde.decode())
    ## --------------------------------------------------------------------
    ## --dim
    print(f'** Dimension file: {args.dim}')
    fn_dim=os.path.join(args.dst,'src/vpsc_src/parameters.dim')
    if args.dim!='src/vpsc_src/parameters.dim':
        print(f'*** args.dim: {args.dim}')
        shutil.copy(args.dim,fn_dim)

    ## --------------------------------------------------------------------
    ## --nphs
    print(f'** Switch <nphmx> in the dimension: given nphs is {args.nphs}')
    chg_kw_in_dim(fn_dim,'nphmx',args.nphs)
    ## --ngr
    ## change ngr
    print(f'** Switch <ngrmx> in the dimension: given ngrs is {args.ngrs}')
    if args.ngrs!=[]:
        ngrs=0
        for i in range(len(args.ngrs)):
            ngrs=ngrs+int(args.ngrs[i])
        chg_kw_in_dim(fn_dim,'ngrmx',ngrs)
    ## --------------------------------------------------------------------
    ## --ielast
    ## change the behavior of the polycrystal
    ## (either elasto-visco-plastic or elastic)
    fn_main=os.path.join(args.dst,'src','vpsc_src','evpsc_main.f')
    if args.ielast:
        str0='ielast=0'
        str1='ielast=1'
        print('\n\n** suppressing plasticity by increasing CRSS **\n\n')
    else:
        str0='ielast=1'
        str1='ielast=0'
    with open(fn_main,'r') as fo:
        oldstring=fo.read()
    newstring=re.sub(str0,str1,oldstring)
    with open(fn_main,'w') as fo:
        fo.write(newstring)

    ## --------------------------------------------------------------------
    ## --ioptfs
    ## Change the ioptfs flag in the main file
    fn_main=os.path.join(args.dst,'src','vpsc_src','evpsc_main.f')
    with open(fn_main,'r') as fo:
        oldstring=fo.read()
    newstring, n = re.subn(
        r"^.*parameter\(ioptfs=.*$",
        f"      parameter(ioptfs={args.ioptfs})",
        oldstring,count=1,flags=re.MULTILINE,)
    if n != 1:
        print("Could not find the ioptfs parameter")
        sys.exit(-1)
        raise ValueError("Could not find the ioptfs parameter")

    with open(fn_main,'w') as fo:
        fo.write(newstring)

    ## --------------------------------------------------------------------
    ## --------------------------------------------------------------------
    ## --iodf
    ## change iodf flag and dim file
    print(f'** Change iodf={args.iodf}')
    chg_iodf(iodf=args.iodf,filename=fn_main)


    ## --------------------------------------------------------------------
    ## --irgvb
    ## change for RGVB (irgvb true or false)
    if args.irgvb:
        str0,str1='irgvb=.false.','irgvb=.true.'
    elif not(args.irgvb):
        str0,str1='irgvb=.true.','irgvb=.false.'

    with open(os.path.join(args.dst,'src','vpsc_src','statev_io.f'),'r') as fo:
        oldstring=fo.read()
    newstring=re.sub(str0,str1,oldstring)
    with open(os.path.join(args.dst,'src','vpsc_src','statev_io.f'),'w') as fo:
        fo.write(newstring)

    ## --------------------------------------------------------------------
    ## Build
    print('** Final dim used in build.py')
    with open(fn_dim,'r') as fo:
        print(fo.read())
    print('** Build')
    print('*** clean test')
    pc=subprocess.Popen(['make','-i','clean'],cwd=os.path.join(args.dst,'umat/src/evpsc_incr'),
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdo,stde=pc.communicate()

    try:
        print(stdo.decode(codec))
        print(stde.decode(codec))
    except:
        print(stdo.decode())

    ## --------------------------------------------------------------------
    ## --idist (distribution version or not)
    print('*** Build test')
    cmd_compile_test='make -i -j test -f Makefile FC=gfortran FFLAGS=-I../../../src/vpsc_src/'
    sys.stdout.write(f'cmd: {cmd_compile_test}')

    pb=subprocess.Popen(cmd_compile_test.split(),
                        cwd=os.path.join(args.dst,'umat/src/evpsc_incr'),
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdo,stde=pb.communicate()
    print('stdo and stde from building test')
    print(stdo.decode())
    print(stde.decode())
    try:
        print('*** try printing make stdo and stde')
        print(stdo.decode(codec))
        print(stde.decode(codec))
    except:
        print('*** make stdo and stde failed')
        pass

    print('*** Run test')
    print(f'*** sys.platform: "{sys.platform}"')
    if sys.platform=='win32':
        pc=subprocess.Popen(['test.exe'],cwd=os.path.join(args.dst,'umat/src/evpsc_incr'),
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    else:
        pc=subprocess.Popen(['./test'],cwd=os.path.join(args.dst,'umat/src/evpsc_incr'),
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    print('*** End of test')
    # pc.wait()
    stdo,stde=pc.communicate()
    print('\n\n*stdo:\n',stdo.decode(codec))
    print('\n\n*stde:\n',stde.decode(codec))


    print('** Compile dEVPSC-SA')
    cmd_compile=['make','-j','evpsc']
    if args.idebug:
        cmd_compile.append('-e')
        cmd_compile.append('ENV=debug')

    pbd=subprocess.Popen(cmd_compile,cwd=args.dst,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdo,stde=pbd.communicate()
    try:print('\n\n*stdo:\n',stdo.decode(codec))
    except:pass
    try:print('\n\n*stde after "make -j evpsc":\n',stde.decode(codec))
    except:pass

    print('** End of build')
    sys.exit(0)
