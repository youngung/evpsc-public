"""
Build (vumat or umat) with Python command environment mostly for debugging.
The plan is to make this module an exclusive method to build UMAT/VUMAT of EVPSC code.

Youngung Jeong
yjeong@changown.ac.kr

Build consists of three steps:
   1. Build prep
   2. Stand-alone build
   3. FE build (UMAT)
"""
import re

def chg_parameter_val(fn=None,parameter=None,val=None,iopt=0):
    """
    Arguments
    ---------
    fn
    parameter
    val
    iopt=0 (0 for value, 1 for string)
    """
    # change parameters
    if iopt==0: new=f'      parameter({parameter}={val})\n'
    if iopt==1: new=f'      parameter({parameter}="{val}")\n'

    with open(fn,'r') as fo:
        lines=fo.readlines()
        for i, line in enumerate(lines):
            line = line.replace(" ","")
            line.lower()
            if f'parameter({parameter}=' in line:
                lines[i]=new
                # print('found an incident')
                break # Applies only to the first incident?

    with open(fn,'w') as fo:
        for i, line in enumerate(lines):
            fo.write(line)

def build_prep(dst,conds,ioptfs,
               idebug=False,irgvb=False,idd=False,
               iregime=2):
    """
    Preparation the location (dst) ready for building EVPSC (either SA or FE).

    Arguments
    ---------
    dst
    conds  : dictionary of parameters specified in the dEVPSC-FE input file.
    ioptfs : 0: IJP 2020, 1: Still affine but with a secant compliance tensor; Applicable only for EVPSC
    idebug : Boolean flag to tell if debug mode is used for compilation
    irgvb  : Boolean flag to tell if RGVB hardening is used (or not).
    idd    : Boolean flag to tell if dis. density-based hardening is used (or not)
    iregime : integer to tell the modeling regime:  0: VPSC
                                                    1: dEVPSC (delta; incremental formulation)
    """
    import os
    import subprocess
    import shutil
    import sys
    iodf=conds['iodf']
    isolver=conds['isolver']
    kwargs=dict(stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    dst_ori=os.getcwd()
    dst_base=dst[::]
    for i in range(3):
        dst_base=os.path.split(dst_base)[0]

    print("*** Build preparations")
    #    * 2-2. Generate dimension files based on the given single crystal file & texture
    cmd=f"python -u pyscripts/gen_dim.py --iregime {iregime} "
    cmd=f"{cmd} --ngrmx {conds['ntotgr']} --nphmx {conds['nph']} "
    cmd=f"{cmd} --fnsx "
    fnsx_all=''
    for iph,_fnsx_ in enumerate(conds['sx']): ## the single crystal file names has been changed.
        __fnsx__=f'sx_{iph}'
        fnsx_target=os.path.join(dst,__fnsx__)
        fnsx_all=f'{fnsx_all} {fnsx_target}'
        try:
            shutil.copy(_fnsx_,fnsx_target)
        except shutil.SameFileError:
            sys.stderr.write( '** Skipped: source and destination are the exact same file')
            sys.stderr.write(f'    source: {_fnsx_} and destination: {fnsx_target}')
        print(f'iph: {iph}, fnsx: {_fnsx_} -> fnsx_target: {fnsx_target}')
        if not os.path.isfile(fnsx_target):
            print(f' os.path.isfile(fnsx_target): {os.path.isfile(fnsx_target)}')
            print(f'** Error: Could not find {fnsx_target}')
            sys.exit(-1)

    cmd=f"{cmd} {fnsx_all} -o src/vpsc_src/parameters.dim"
    print(f'*** <gen_dim> command used for build: {cmd}')
    p=subprocess.Popen(cmd.split(),cwd=dst_base,**kwargs)

    stdo,stde=p.communicate()
    if p.poll()!=0:
        print(stdo.decode())
        print(stde.decode())
        print(f'** Error: cmd <{cmd}> failed')
        print('**** dst     :',dst)
        print('**** dst_base:',dst_base)
        sys.exit(0)

    os.chdir(dst_base) # I am not sure if this is necessary.

    ## irgvb flag
    print('**** Build 1: RGVB or DD hardening flag')
    print('***** 1-1: IRGVB flag (RGVB hardening model or not)')
    with open(os.path.join(dst_base,'src','vpsc_src','statev_io.f'),'r') as fo:oldstring=fo.read()
    if irgvb:
        newstring=re.sub('irgvb=.false.','irgvb=.true.',oldstring)
    else:
        newstring=re.sub('irgvb=.true.','irgvb=.false.',oldstring)
    with open(os.path.join(dst_base,'src','vpsc_src','statev_io.f'),'w') as fo:fo.write(newstring)
    print('***** 1-2: IDD flag (DD hardening model or not)')
    with open(os.path.join(dst_base,'src','vpsc_src','statev_io.f'),'r') as fo:oldstring=fo.read()
    if idd:
        newstring=re.sub('idd=.false.','idd=.true.',oldstring)
    else:
        newstring=re.sub('idd=.true.','idd=.false.',oldstring)
    with open(os.path.join(dst_base,'src','vpsc_src','statev_io.f'),'w') as fo:fo.write(newstring)

    ## iodf change implicit or explicit FE solver build
    print('**** Build 2: iodf change in <sdvini.f> for implicit FE solver or in <vumat.f> for explicit FE solver')
    if isolver==0: ## implicit solver (UMAT)
        ## sdvini.f
        chg_parameter_val(fn=os.path.join(dst,'sdvini.f'),parameter='iodf',val=iodf,iopt=0)
    elif isolver==1: ## explicit solver (VUMAT)
        ## inside <vumat.f> from which <evpsc_input_fe> is called.
        chg_parameter_val(fn=os.path.join(dst,'vumat.f'),parameter='iodf',val=iodf,iopt=0)
    chg_parameter_val(fn=os.path.join(dst,'..','..','..','src','vpsc_src','evpsc_main.f'),
                      parameter='iodf',val=iodf,iopt=0)

    ## ioptfs
    chg_parameter_val(fn=os.path.join(dst,'..','..','..','src','vpsc_src','evpsc_main.f'),
                      parameter='ioptfs',val=ioptfs,iopt=0)
    if isolver==0:
        chg_parameter_val(fn=os.path.join(dst,'umat.f'),
                          parameter='ioptfs',val=ioptfs,iopt=0)
    elif isolver==1:
        chg_parameter_val(fn=os.path.join(dst,'vumat.f'),
                          parameter='ioptfs',val=ioptfs,iopt=0)

    ## environment
    fn=os.path.join(dst_ori,'env','path_to_env.txt')
    with open(fn,'r') as fo:
        path_to_env=fo.read().split('\n')[0]

    if idebug:
        fn_ref=os.path.join(dst_ori,'env','win86_64_debug.env')
    else:
        fn_ref=os.path.join(dst_ori,'env','win86_64_yj.env')

    print('idebug:',idebug)
    print(f'fn_ref: {fn_ref}')
    print(f'dst_ori: {dst_ori}')
    print(f'path_to_env: {path_to_env}')
    shutil.copy(fn_ref,path_to_env)


def build_sa(dst):
    """
    Build stand-alone code.

    Arguments
    ---------
    dst    : location in which build is performed.
    """
    import subprocess, sys
    kwargs=dict(stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=dst)

    ## stand-alone build
    print('**** Build 4: Stand-alone build if necessary')
    p_sa=subprocess.Popen('make -i clean'.split(),**kwargs)
    stdo,stde=p_sa.communicate()

    if sys.platform=='win32':
        cmd_sa_build='make -j evpsc -f win_ifort.mk'
    elif sys.platform=='linux' or sys.platform=='linux2':
        cmd_sa_build='make -j evpsc -f nix_ifort.mk'
    elif sys.platform=='darwin':
        cmd_sa_build='make -j evpsc -f mac_ifort.mk'
    else:
        raise IOError(f'** Error: Unexpected platform {sys.platform}')

    p_sa=subprocess.Popen(cmd_sa_build.split(),**kwargs)
    stdo,stde=p_sa.communicate()
    print(f'***** cmd_sa_build: {cmd_sa_build}')

    try:
        print(stdo.decode(errors='ignore'))
        print(stde.decode(errors='ignore'))
    except:
        decoder1=chardet.detect(stdo)
        decoder2=chardet.detect(stde)
        print(f'decoder1: {decoder1}')
        print(f'decoder2: {decoder2}')

def build_fe(dst,isolver):
    """
    Build UMAT or VUMAT (as well as the stand-alone version if isa is True).
    Building VUMAT is for testing only. Actual running of VUMAT requires
    building it on the fly by passing the source code to Abaqus.

    Arguments
    ---------
    dst
    isolver
    """
    import sys,os
    import subprocess
    kwargs=dict(stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=dst)
    print('**** Build 5: make test')
    print('*****      5-1: make clean')
    p=subprocess.Popen(['make','-i','clean'],**kwargs)
    stdo,stde=p.communicate()
    stdo,stde=p.communicate()
    print(f'p.poll(): {p.poll()}')

    print('*****      5-1: make test')
    p=subprocess.Popen(['make','-i','test','FC=ifort'],**kwargs)
    stdo,stde=p.communicate()
    print(f'p.poll(): {p.poll()}')
    try:
        print('stdo')
        print(stdo.decode(errors='ignore'))
        print('stde')
        print(stde.decode(errors='ignore'))
    except:
        decoder1=chardet.detect(stdo)
        decoder2=chardet.detect(stde)
        print(f'decoder1: {decoder1}')
        print(f'decoder2: {decoder2}')

    print('**** Build 6: run test')
    if sys.platform=='win32':
        if isolver==0:
            p=subprocess.Popen('test.exe',shell=True,**kwargs)
        elif isolver==1:
            p=subprocess.Popen(['test.exe','1','1'],**kwargs)
    else:
        if isolver==0:
            p=subprocess.Popen('./test',shell=True,**kwargs)
        elif isolver==1:
            p=subprocess.Popen(['./test','1','1'],**kwargs)

    stdo,stde=p.communicate()
    p.wait()

    try:
        print('stdo')
        print(stdo.decode(errors='ignore'))
        print('stde')
        print(stde.decode(errors='ignore'))
    except:
        decoder1=chardet.detect(stdo)
        decoder2=chardet.detect(stde)
        print(f'decoder1: {decoder1}')
        print(f'decoder2: {decoder2}')

    print(f'p.poll(): {p.poll()}')
    if p.poll()!=0:
        print(f'sys.platform: {sys.platform}')
        print(f'solver      : {isolver}')
        print(f'dst         : {dst}')
        sys.stderr.write('** Error! test resulted in non-zero results')
        sys.exit(-1)

    print('**** Build 7: prebuild')
    fnlog=os.path.join(os.getcwd(),'log.txt')
    print(f"log: {fnlog}")

    with open(fnlog,'w') as fo:
        sys.stdout.write('** diffraction build\n')
        fo.write('** diffraction build\n')
        p_=subprocess.Popen(['make','-i','diffraction','FC=gfortran','FFLAGS=-I../../../src/vpsc_src/'],**kwargs)
        stdo_,stde_=p_.communicate()

        if isolver==0:
            p=subprocess.Popen(
                ['python','-u','prebuild_vumat.py','-i','2'],
                **kwargs)
        elif isolver==1:
            p=subprocess.Popen(
                ['python','-u','prebuild_vumat.py','-i','3'],
                **kwargs)
        stdo,stde=p.communicate()
        print('print stdo & stde')
        fo.write(stdo.decode(errors='ignore'))
        fo.write(stde.decode(errors='ignore'))
        try:
            print(stdo.decode(errors='ignore'))
            print(stde.decode(errors='ignore'))
        except:
            decoder1=chardet.detect(stdo)
            decoder2=chardet.detect(stde)
            print(f'decoder1: {decoder1}')
            print(f'decoder2: {decoder2}')


        print('-----------------------------------------')
        print(f'sys.platform: {sys.platform}')

        if isolver==0:
            print('** Create UMAT in <run.py>')
            if sys.platform=='win32':

                p=subprocess.Popen(['make','-i'],**kwargs)
                print('** diffraction build')
                stdo,stde=p.communicate()

                stdo_msg=stdo.decode(errors='ignore')
                stde_msg=stde.decode(errors='ignore')

                error_message='make: [umat-std.obj] Error 1 (ignored)'
                if error_message in stde_msg:
                    print(stdo.decode(errors='ignore'))
                    print(stde.decode(errors='ignore'))

                    sys.stderr.write('** Error occured during build umat-std.obj\n')
                    return -1


            elif sys.platform=='linux' or sys.platform=='linux2':
                p=subprocess.Popen(['make','-i','umat-std.o'],**kwargs)
                stdo,stde=p.communicate()
                try:
                    print(stdo_.decode(errors='ignore'))
                    print(stde_.decode(errors='ignore'))
                except:
                    decoder1=chardet.detect(stdo)
                    decoder2=chardet.detect(stde)
                    print(f'decoder1: {decoder1}')
                    print(f'decoder2: {decoder2}')

            elif sys.platform=='darwin':
                print("** Cannot compile umat-std.o in Mac")
                sys.exit(-1)
            else:
                print('** Unexpected sys.platform')
                sys.exit(-1)


        elif isolver==1:
            print('** VUMAT is not precompiled to an object,')
            print('   but will be compiled only on the fly')
            print('   when submitting the actual abaqus job')
            print('-----------------------------------------')

            if sys.platform=='win32':
                p=subprocess.Popen(['make','-i','vumat-xplD.obj'],**kwargs)
            elif sys.platform=='linux' or sys.platform=='linux2':
                p=subprocess.Popen(['make','-i','vumat-xpl.o'],**kwargs)
            elif sys.platform=='darwin':
                print("** Cannot compile vumat in Mac")
                sys.exit(-1)
            else:
                print('** Unexpected sys.platform')
                sys.exit(-1)

            stdo,stde=p.communicate()
        print('\nprinting stdo&stde')

        try:
            print(stdo.decode(errors='ignore'))
            fo.write(stdo.decode(errors='ignore'))
            print(stde.decode(errors='ignore'))
            fo.write(stde.decode(errors='ignore'))
        except:
            decoder1=chardet.detect(stdo)
            decoder2=chardet.detect(stde)
            print(f'decoder1: {decoder1}')
            print(f'decoder2: {decoder2}')

    print('**** <run.build> completed')
    return 0

import sys
if __name__=='__main__':
    sys.stdout.write("""* beginning <build_umat.py>\n""")
    import os
    import numpy as np
    import argparse

    import tempfile
    import shutil
    import subprocess

    sys.stdout.write("""--------------------------------------------------------------------------\n""")
    sys.stdout.write("""       Build umat or vumat of dEVPSC-FE code for Abaqus simulation\n""")
    sys.stdout.write("""          Youngung Jeong; yjeong@changwon.ac.kr\n""")
    sys.stdout.write("""--------------------------------------------------------------------------\n\n""")
    ## 1. Argument parsing
    parser = argparse.ArgumentParser()


    ## Common arguments shared with <run.py>
    # parser.add_argument('--fn',type=str,help='Input file to <run.py>; examples are relocated in <umat/src/evpsc_incr/tests> folder. The running conditions given in this file is used but can be overriden by arguments passed.',default='umat.in')
    parser.add_argument('--ngrs',nargs='+',default=[0],help='The total number of grains per each phase within each FE integration point.')
    parser.add_argument('--iodf',type=int,help='0: discrete; 1: odf (mtex), 2: random, 3: EBSD, 4: Euichan, 5: shuffle')
    parser.add_argument('--idd', action='store_true',help='Boolean flag to switch on the DD hardening model; switch in <statev_io.f>')
    parser.add_argument('--irgvb', action='store_true',help='Boolean flag to switch on the RGVB-DD hardening model; switch in <statev_io.f>')
    parser.add_argument('--idebug', action='store_true',help='Boolean flag to switch on debugging options')
    parser.add_argument('--iregime',type=int,default=-1,help='1: VPSC; 2: delta-EVPSC (incremental formulation)')
    parser.add_argument('--ioptfs',type=int,default=-1,help='0: IJP 2020, 1: Still affine but with a secant compliance tensor; Applicable only for EVPSC')

    parser.add_argument('--isolver',default=0,type=int,help='abaqus/standard (0; implicit) or abaqus/explicit (1; explicit)')

    parser.add_argument('--pathwork',type=str,default=os.getcwd(),help='Name of path, in which the main source code either <umat.f> or <vumat.f> is located.')

    ## Arguments that are read from <umat input file> for the case of <run.py>
    parser.add_argument('--fnsx',nargs='+',default=['../../../examples/bench_MMM2018/sts_1e-3.sx'],help='Single crystal file for each phase')
    parser.add_argument('--idiff',nargs='+',default=[0],help='IDIFF; The flag to do diffraction simulation for each phase')
    parser.add_argument('--diff',nargs='+',default=[''],help='diff condition file for each phase')
    parser.add_argument('--wph',nargs='+',default=[1.0],help='weight fraction for each phase')
    parser.add_argument('--nph',default=1,type=int,help='The total number of phases')
    parser.add_argument('--irollback',action='store_true',help='Whether or not to roll-back changes introduced over the course of build')

    args=parser.parse_args()

    ## testing
    if args.iregime==-1:
        sys.stderr.write('** iregime should be specified')
        sys.exit(-1)
    if args.iregime!=1 and args.iregime!=2 and args.iregime!=3:
        sys.stderr.write('** iregime should be among 1, 2, and 3')
        sys.exit(-1)

    if args.iregime==2 and args.ioptfs==-1:
        sys.stderr.write('** For EVPSC, (iregime==1), ioptfs should be specified')
        sys.exit(-1)


    # Directory that is three level below than given pathwork.
    dst_base=args.pathwork[::]
    for i in range(3):
        dst_base=os.path.split(dst_base)[0]

    conds=dict()

    conds.update(isolver=args.isolver,
               text=None,
               iodf=args.iodf,
               sx=args.fnsx,
               idiff=args.idiff,
               diff=args.diff,
               abqinp=None,
               probs=None,
               ipx=0,
               pxfn=None,
               nph=args.nph,
               wph=args.wph)

    ## update total number of discrete grains.
    args.ngrs=np.array(args.ngrs,dtype='int')
    ntotgr=args.ngrs.sum()
    conds.update(ntotgr=ntotgr)

    print('*******************************')
    print('** copy sx to sx_iph file')
    fnsx_news=[]
    for iph,_fnsx_ in enumerate(conds['sx']): ## the single crystal file names has been changed.
        fnsx_new=f'sx_{iph}'
        fnsx_old=f'{os.path.join(args.pathwork,_fnsx_)}'
        print(f'    copying {fnsx_old} to {fnsx_new}')
        shutil.copy(fnsx_old,fnsx_new)
        fnsx_news.append(fnsx_new)

    conds.update(sx=fnsx_news)
    print('*******************************')
    print('** Build begins')
    print('** Build prep')
    build_prep(dst=args.pathwork,conds=conds,
               ioptfs=args.ioptfs,
               idebug=args.idebug,irgvb=args.irgvb, idd=False,
               iregime=args.iregime)

    # print('** Build stand-alone')
    # build_sa(dst=dst_base)
    print('** Build FE')
    isig=build_fe(dst=args.pathwork,isolver=args.isolver)
    print('** Build prep ends')
    print('*******************************')

    if args.irollback:
        ## roll-back changes introduced by prebuild_vumat.py
        print('\n\n** Roll-back changes introduced by <prebuild_vumat.py>')
        cmd='python -u prebuild_vumat.py'
        if args.isolver==0:  cmd=f'{cmd} -i 4'
        if args.isolver==1:  cmd=f'{cmd} -i 5'
        print(f'** cmd for Roll-back changes {cmd}')
        p=subprocess.Popen(
            cmd.split(),cwd=args.pathwork,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdo,stde=p.communicate()
        print('\n** stdo and stde after roll-back')
        print(stdo.decode(errors='ignore'))
        print(stde.decode(errors='ignore'))
        print('-----------------------------------')

    if isig==0:
        print('** Build based on <build_umat.py> has been completed')
    else:
        print('** Build based on <build_umat.py> failed')
