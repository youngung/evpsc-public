"""
build test inside umat folder and run the test.
"""
import os
import subprocess
def run_popen(cmd,cwd=os.getcwd(),echo=False):
    """
    """
    import sys
    kws=dict(stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p=subprocess.Popen(cmd,cwd=cwd,**kws)

    stdo,stde=p.communicate()
    if echo:
        print('\n** Echo begins')
        print(stdo.decode())
        print(stde.decode())
        print(f'poll: {p.poll()}')
        print('** Echo ends\n')
    return p.poll()

if __name__=='__main__':
    import argparse
    import sys
    import tempfile
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--iopt',type=int,
        help='0: compiles with gfortran; 1: compiles with ifort',default=0)
    args=parser.parse_args()

    ## check if the compiler available.
    if args.iopt==0:
        FC='gfortran'
        p=subprocess.Popen([FC,'--version'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    elif args.iopt==1:
        FC='ifort'
        p=subprocess.Popen([FC,'/help'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    stdo,stde=p.communicate()
    # print(stdo.decode()); print(stde.decode())
    sys.stdout.write(f'** Is {FC} available? {p.poll()}')
    if p.poll()!=0:
        sys.stderr.write(f'** Error: could not find {FC}')
        sys.exit(-1)

    ## Clean test
    fnmake='Makefile' # 'makefile_gfortran.mk'

    clean_cmd=f'make -i clean -f {fnmake}'
    run_popen(clean_cmd.split(),cwd='umat/src/evpsc_incr',echo=False)

    compile_cmd=f'make -i test -f {fnmake} FC={FC}'
    if args.iopt==0:
        compile_cmd=f'{compile_cmd} FFLAGS=-I../../../src/vpsc_src/'
    run_popen(compile_cmd.split(),cwd='umat/src/evpsc_incr',echo=True)

    ## -- check system.platform
    sys.stdout.write(f'** sys.platform: {sys.platform}')
    if sys.platform=='win32':
        sys.stdout.write('\t A windows system detected')
        run_cmd='./test.exe'
    else: ## Darwin (Mac) or Linux
        sys.stdout.write('\t A non-windows system detected')
        run_cmd='./test'


    sys.stdout.write(f'\n** run_cmd: {run_cmd}')
    if False:
        poll = run_popen(cmd=run_cmd,cwd='umat/src/evpsc_incr',echo=True)
    else:
        home=os.getcwd()
        os.chdir('umat/src/evpsc_incr')
        p=subprocess.Popen(run_cmd,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        stdo,stde=p.communicate()
        os.chdir(home)
        sys.stdout.write(f'\n** after {run_cmd}')
        sys.stdout.write('\n** stdo\n')
        sys.stdout.write(stdo.decode())
        sys.stderr.write('\n** stde\n')
        sys.stderr.write(stde.decode())
        sys.stdout.write(f'** p.poll : {p.poll()}')
