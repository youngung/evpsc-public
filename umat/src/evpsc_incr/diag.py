## This script is to test <input> files used for UMAT in conjunction with <test_nstatv_depvar.f>
## Also can change specific lines of all input files by giving some tweaks.

from glob import glob
import __future__
import re
import time


def read_input(fn):
    """
    Arguments
    ---------
    fn: Abaqus Input file <*.inp>

    Returns
    -------
    Find
    """
    with open(fn,'r') as fo:
        cont=fo.read()
        cont=re.split('Depvar',cont,flags=re.IGNORECASE)[1]
        lines=cont.split('\n')

    return int(lines[1].split(',')[0])

def check_sdvini(fn):
    with open(fn,'r') as fo:
        cnt=fo.read()
    if "*initial conditions, type=solution, user" in cnt:
        return True
    else:
        return False

def main(q=False):
    """
    q = Quiet
    """
    #print('')
    _fns_=glob('*.inp')
    dat=[]
    fns_regular_abq_inp=[]
    fns_restart=[]
    for i in range(len(_fns_)):

        try:
            ## Files with DEPVAR included.
            n=read_input(_fns_[i])
            dat.append(n)
            fns_regular_abq_inp.append(_fns_[i])
        except:
            ## see if the file is 'restart' input.
            try:
                fn=_fns_[i]
                # print('try again %s'%fn)
                with open(fn,'r') as fo:
                    cont=fo.read()
                    cont=re.split(r'Restart, read',cont,flags=re.IGNORECASE)[1]
                    lines=cont.split('\n')
                # print('** Restart file found')
                fns_restart.append(fn)
            except:
                # print(f'failed to test {_fns_[i]}') # Python 3
                # if not(q): print('Failed to test <%s>'%_fns_[i]) # Python 2 with Future or Python 3
                pass



    fn='list_abq_inputs.txt'
    with open(fn,'w') as fo:
        fo.write('Updated date&time:%s\n'%time.asctime()) # Python 2 with Future or Python 3
        fo.write(' %i  Total # of tested input files \n'%len(dat)) # python 2 with Future or Python 3
        for i in range(len(fns_regular_abq_inp)):
            fo.write('%s %s\n'%(dat[i],fns_regular_abq_inp[i])) # Python 2 with Future or Python 3
        if len(fns_restart)>0:
            fo.write('-- below is the list of restart input files\n')
            for i in range(len(fns_restart)):
                fo.write('%s\n'%(fns_restart[i])) # Python 2 with Future or Python 3

    # if not(q):
    #     #print('')
    #     with open(fn, 'r') as fo:
    #         print(fo.read())

    for i in range(len(_fns_)):
        if not(check_sdvini(_fns_[i])):
            print('Warning! <initial condition> is not found')

    return fns_regular_abq_inp,fns_restart

def fix_depvar(fns,nvar):
    """
    Change Depvar in Abaqus input files.
    """
    for i in range(len(fns)):
        # print(fns[i])
        with open(fns[i],'r') as fo:
            orig_cont=fo.read() ## original contents

        part_head=re.split('Depvar',orig_cont,flags=re.IGNORECASE)[0]
        part_tail=re.split('Depvar',orig_cont,flags=re.IGNORECASE)[1]
        tail_lines=part_tail.split('\n')[2:]
        part_new_tail='Depvar\n%i'%(nvar)
        new_cont='%s%s'%(part_head,part_new_tail)
        for j in range(len(tail_lines)):
            new_cont='%s\n%s'%(new_cont,tail_lines[j])

        with open(fns[i],'w') as fo:
            fo.write(new_cont)


def fix_elem_output(fns,keys):
    """
    Change element output in Abaqus input files.

    Arguments
    ---------
    fns
    keys
    """
    statement_ref='*Element Output'
    nd=len(statement_ref)

    for i in range(len(fns)):
        with open(fns[i],'r') as fo:
            orig_cont=fo.read() ## original contents
            orig_lines=orig_cont.split('\n')

        if statement_ref in orig_cont:
            nid=[]
            for j in range(len(orig_lines)):
                l=orig_lines[j]
                if l[:nd]==statement_ref:
                    nid.append(j+1)

            for j in range(len(nid)):
                dd=''
                for k in range(len(keys)):
                    dd='%s%s,'%(dd,keys[k].strip())
                orig_lines[nid[j]]=dd

            with open(fns[i],'w') as fo:
                for j in range(len(orig_lines)):
                    fo.write('%s\n'%orig_lines[j])


def fix_node_output(fns,keys):
    """
    Change nodal output in Abaqus input files.

    Arguments
    ---------
    fns
    keys
    """
    statement_ref='*Node Output'
    nd=len(statement_ref)

    for i in range(len(fns)):
        # print(fns[i])
        with open(fns[i],'r') as fo:
            orig_cont=fo.read() ## original contents
            orig_lines=orig_cont.split('\n')

        if statement_ref in orig_cont:
            nid=[]
            for j in range(len(orig_lines)):
                l=orig_lines[j]
                if l[:nd]==statement_ref:
                    nid.append(j+1)

            for j in range(len(nid)):
                dd=''
                for k in range(len(keys)):
                    dd='%s%s,'%(dd,keys[k].strip())
                orig_lines[nid[j]]=dd

            with open(fns[i],'w') as fo:
                for j in range(len(orig_lines)):
                    fo.write('%s\n'%orig_lines[j])


def set_mass_scale(fns,msf):
    """
    insert mass scale options for steps
    """
    statement_ref="*Dynamic, Explicit"
    nd=len(statement_ref)

    for i in range(len(fns)):
        with open(fns[i],'r') as fo:
            orig_cont=fo.read() ## original contents
            orig_lines=orig_cont.split('\n')

        if statement_ref in orig_cont:
            nid=[]
            for j in range(len(orig_lines)):
                l=orig_lines[j]
                if l[:nd]==statement_ref:
                    nid.append(j+2)

            statement=f"*FIXED MASS SCALING, FACTOR={msf}"
            # statement="*FIXED MASS SCALING, TYPE=UNIFORM, DT=1e-4"
            #statement="*FIXED MASS SCALING, FACTOR=1e2\n*VARIABLE MASS SCALING"

            for j in range(len(nid)):
                orig_lines.insert(nid[j]+j,statement)

            ## un comment-out belo to actually rewrite
            """
            for j in range(len(nid)):
                dd=''
                for k in range(len(keys)):
                    dd='%s%s,'%(dd,keys[k].strip())
                orig_lines[nid[j]]=dd



            """
            with open(fns[i],'w') as fo:
                for j in range(len(orig_lines)):
                    fo.write('%s\n'%orig_lines[j])



if __name__=='__main__':
    """
    Check depvar keywords in Abaqus input files
    (And optionally) change *depvar keywords therein

    python diag.py -n <nstatv> -i [0,1] -eok <list of eok> -nok <list of nok> -msf <mass scale factor>
    if <nstatv> is not given, it does not update.
    """
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',type=int,help='standard (0; UMAT) or explicit (1; vumat)',default=0)
    parser.add_argument('-n',type=int,help='Number of state variables',default=0)
    parser.add_argument('-msf',type=float,help='mass scale factor',default=1)

    args = parser.parse_args()

    if args.i==0:
        default_eok=['COORD','LE','S','UVARM','ER','NFORC']
        default_nok=['U','UR','COORD','RF','RT','RM','TF']
    elif args.i==1:
        default_eok=['COORD','LE','S','ER','NFORC']
        default_nok=['U','UR','COORD','RF','RT','RM']

    parser.add_argument('-eok',nargs='+',help='Element output key list',
                        default=default_eok)
    parser.add_argument('-nok',nargs='+',help='Nodal output key list',
                        default=default_nok)

    args = parser.parse_args()

    fns,fns_restart=main(q=False) ## regular inputs

    #
    # fns=['3pb_last.inp']
    # fns=[]

    # print('From now on, only 3pb_last.inp is modified by "test.exe"')

    if args.n!=0:
        fix_depvar(fns,args.n) ## regular inputs
        main(q=True)

        # This check-up before writing 'umat-nstv.f' is to
        # avoid redundant build by Makefile, which detects the
        # necessity of rebuilding 'umat-std.obj' based on
        # the time of modification (or creation) made in dependent files.
        # 'umat-nstv.f' is one of such dependents.

        nstv=-1 ##
        if os.path.isfile('umat-nstv.f'):
            with open('umat-nstv.f','r') as fo:
                l=fo.readlines()[1]
                nstv=int(l.split('=')[-1].split(')')[0])

        if nstv!=args.n:
            with open('umat-nstv.f','w') as fo:
                fo.write('      integer nstatv\n')
                fo.write('      parameter(nstatv=%i)'%args.n)








    #print('eok:',args.eok)
    #print('nok:',args.nok)

    #for i in range(len(fns)):
    fix_elem_output(fns,args.eok)
    fix_node_output(fns,args.nok)

    if len(fns_restart)>0:
        fix_elem_output(fns_restart,args.eok)
        fix_node_output(fns_restart,args.nok)



    ## scale mass
    if args.i==1:
        # print('diag.py runs for abaqus/explicit solver')
        # print('mass scale factor:',args.msf)
        set_mass_scale(fns,args.msf)

    ##
    if args.i==0:
        print('** from UMAT_v2.2.0 onward, one should use <initial conditions> keyword for Abaqus/standard')
        _fns_=glob('*.inp')
        for i in range(len(_fns_)):
            if not(check_sdvini(_fns_[i])):
                print('** <initial condtion> could not be found')
                #raise IOError('input file %s does not use <initial condition>'%_fns_[i])
