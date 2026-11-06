"""
- This script, <run.py>, is to help running dEVPSC-FE simulation easily.
- This script is used to build the necessary objective file for Abaqus
- This script helps to modify the Abaqus input file.
"""
import os, sys
import datetime
import shutil
import tarfile
import chardet
import numpy as np
import re
from build_umat import chg_parameter_val


__umat_in_template__="""** input to run dEVPSC-FE run via UMAT
* iodf (0: discrete ori; 1: ODF in 5deg grid (mtex format); 2: iso; 3: EBSD (ang fom); 4: spread of an ideal component
2
* nph (number of phases)
2
* texture/ODF file;; eddq_push/EDDQ100
../../../matData/vpscData/EDDQ/eddq.odf
../../../matData/vpscData/EDDQ/eddq.odf
* sx file
eddq_push/eddq_yj.sx
eddq_push/eddq_yj.sx
* idiff (0: no diffraction; 1: diffraction)
0
* diff file (dummy if idiff = 0 is given in the above line)
matData/lc/bcc.dif
* abaqus input file
inps/UMAT_DEV_VEL_0001.inp
* noel and npt subjected to 'lookup' (could be dummy unless <isa> flag is on)
1
1 1
* ipx
0
* pxfn (dummy if ipx is 0)
dummy
"""

if sys.platform=='win32':
    vers=['2020','2021']
    ABQ=None
    for ver in vers:
        fn=f'C:/SIMULIA/EstProducts/{ver}/win_b64/code/bin/ABQLauncher.exe'
        if os.path.isfile(fn):
            ABQ=fn
    if type(ABQ)==type(None):
        raise IOError('** Error: could not find ABQLauncher.')
        sys.exit(-1)
elif sys.platform=='linux' or sys.platform=='linux2':
    ABQ="/apps/commercial/abaqus/Commands/abq2020hf4"
elif sys.platform=='darwin':
    print("** Warning: You are on Mac and Abaqus is not available")

def read_in(fn='umat.in'):
    """
    Read simulation conditions from given file

    Argument
    --------
    fn = 'umat.in'

    Returns
    -------
    A dictionary that contains keyword-value pairs required to run EVPSC-FE
    """
    import os
    import subprocess
    if not os.path.isfile(fn):
        print(f'** Error could not find {fn}')
        print('\n\n\n--------------------------------------------')
        print('use something like below:')
        print(__umat_in_template__)
        print('--------------------------------------------\n\n\n')
        sys.exit(-1)

    with open(fn,'r') as fo:
        lines=fo.readlines()
    for i, line in enumerate(lines):
        lines[i]=line.split('\n')[0]

    ## implicit (or explicit)
    i=2
    isolver=int(lines[i]) ## 0: UMAT or 1: VUMAT
    # print(f'isolver:{isolver}')

    ## Number of phases
    i+=2
    # print(f'lines: {lines[i]}')
    nph = int(lines[i])
    # if nph!=nphmx:
    #     print('**NPH and NPHMX are not equal to each other')
    #     sys.exit(-1)

    text=[]
    sx=[]
    idiff=[]
    diff=[]

    ## weight fraction
    i+=2
    wph = np.array(lines[i].split(),dtype='float')
    ## use of discrete orientation to generate texture on the fly.
    i+=2
    iodf = int(lines[i].split()[0])
    ## texture
    i+=1
    for iph in range(nph):
        i+=1
        text.append(lines[i])
    ## sx
    i+=1
    for iph in range(nph):
        i+=1
        sx.append(lines[i])
    ## diffraction flag
    i+=1
    for iph in range(nph):
        i+=1
        idiff.append(lines[i])
    ## diffraction
    i+=1
    for iph in range(nph):
        i+=1
        diff.append(lines[i])
    ## abaqus inp
    i+=2
    abqinp = lines[i]

    ## elements and integration points to look up
    i+=2
    nprob=int(lines[i])
    probs=[]
    for j in range(nprob):
        i+=1
        noel,npt=lines[i].split()
        noel=int(noel)
        npt=int(npt)
        probs.append([noel,npt])
    ## ipx
    i+=2
    ipx=int(lines[i])
    ## pxfn
    i+=2
    pxfn=lines[i]

    if ipx==1 or ipx==2:
        if not os.path.isfile(pxfn):
            raise IOError(f"** Could not find the {pxfn}")
        else: pass
        if ipx==2 and True:
            ## figure out if all the files listed in the given pxfn
            ## exists
            with open(pxfn) as fo:
                lines_px=fo.readlines()[1:]
                for fn in lines_px:
                    fn=fn.split('\n')[0]
                    if not(os.path.isfile(fn)):
                        print('** something is wrong')
                        print(f'could not find {fn}')
                        print('fn below:')
                        print(fn)
                        sys.exit(0)

    elif ipx==0:
        pass
    else:
        print('ipx:',ipx)
        print('type(ipx):',type(ipx))
        raise IOError(f"** Could not understand ipx {ipx}")

    return dict(isolver=isolver,
                text=text,
                iodf=iodf,
                sx=sx,
                idiff=idiff,
                diff=diff,
                abqinp=abqinp,
                probs=probs,
                ipx=ipx,
                pxfn=pxfn,
                nph=nph,
                wph=wph)

def check_abqinp(abqinp,isolver):
    ## test if initial conditions is read only for standard solver
    print(f'   * isolver {isolver}')
    if isolver==0:
        with open(abqinp,'r') as fo:
            lines_abq_inp=fo.read().split('\n')
        icheck_init_sol=False
        print(f'   * I am checking initial condition')
        for l in lines_abq_inp:
            if '*initialconditions,type=solution,user' in l.replace(' ','').lower():
                icheck_init_sol=True

        print(f'   * icheck_init_sol: {icheck_init_sol}')

        if not(icheck_init_sol):
            print(f'\n\n*** Failed to find <*initial conditions, type=solution, user> in\n{abqinp}')
            sys.exit(-1)


def prep_vpscin_nph(dst,idiffs,textures,sxs,diffs,ngrs,iodf,
                    ref='matData/vpsc7.in_generic_nurion_dp'):
    """
    arguments
    ---------
    dst
    idiffs
    textures
    sxs
    diffs
    ngrs
    iodf
    ref
    """
    import re
    import os,shutil,subprocess

    ## prepare "vpsc7.in" file
    with open(ref,'r') as fo:
        vpscin=fo.read()
        for iph in range(len(textures)):
            vpscin=vpscin.replace(fr'texture_file_generic_path_{iph}',os.path.join(dst,f'texture_{iph}'))
            vpscin=vpscin.replace(fr'sx_file_generic_path_{iph}',os.path.join(dst,f'sx_{iph}'))
            vpscin=vpscin.replace(fr'diffraction_file_generic_path_{iph}',os.path.join(dst,f'diff_{iph}'))

        ## if diffraction is switched on
        lines=vpscin.split('\n')
        for iph in range(len(textures)):
            lines[14+13*iph]=f'{idiffs[iph]}  ## Diff switched on'

        vpscin=''
        for l in lines:
            vpscin=f'{vpscin}{l}\n'
        vpscin=vpscin[:-1]

    with open(os.path.join(dst,'vpsc7.in'),'w') as fo:
        fo.write(vpscin)

    ##
    print('cwd:',os.getcwd())
    for iph in range(len(textures)):
        shutil.copy(textures[iph],os.path.join(dst,f'texture_{iph}'))
        shutil.copy(sxs[iph],os.path.join(dst,f'sx_{iph}'))

        if idiffs[iph]!=0:
            shutil.copy(diffs[iph],os.path.join(dst,f'diff_{iph}'))

        ## reduce the number of grains if ngrs other than 0 is given

        if ngrs[iph]!=0 and iodf==0:
            if len(textures)>1:
                print('prep_vpscin_nph is not ready for iodf==0 and multiphase case.')
                sys.exit(-1)
            texture_filename=os.path.join(dst,f'texture_{iph}')
            with open(texture_filename,'r') as fo:
                ngr_ori=int(fo.read().split('\n')[3].split()[1])

            fntx_target=f"{os.path.join(dst,'texture')}_{iph}"
            if ngr_ori>ngrs[iph]:
                print( "** Reducing the number of grains in texture")
                print(f"    from {ngr_ori} -> {ngrs[iph]}")
                cmd=f"python -u red_tex.py -i {texture_filename}"
                cmd=f"{cmd} -o {fntx_target}"
                cmd=f"{cmd} -n {ngrs[iph]}"

                print(f"cmd: {cmd}")
                subprocess.run(cmd.split())
            elif ngr_ori<ngrs[iph]:
                print( "** Cannot reducing the number of grains in texture")
                print(f" The original file has fewer orientations than requested")
                print(f" No. of grains in original file:  {ngr_ori}")
                print(f" No. of requested grains:         {ngrs[iph]}")
                sys.exit(-1)
            elif ngr_ori==ngrs[iph]:
                print( "** Cannot reducing the number of grains in texture")
                print(f"  Both original and the requested numbers are: {ngr_ori}")







def archive(hashcode,dst,ngrs):
    """
    Creates a tar file containing the necessary files to run FE simulations,
    and extracts the contents in the target location (dst).

    Arguments
    ---------
    hashcode
    dst
    ngrs
    """
    import subprocess
    from glob import glob
    home=os.getcwd()
    dst_base=dst[::]
    for i in range(3): dst_base=os.path.split(dst_base)[0]

    ## create src/vpsc_src
    os.makedirs(os.path.join(dst_base,'src','vpsc_src'))
    # shutil.copy(dim,os.path.join(dst_base,'src','vpsc_src','parameters.dim'))

    ## create tar
    home_base=home[::]
    for i in range(3):
        home_base=os.path.split(home_base)[0]
    os.chdir(home_base)

    fns=glob('umat/src/evpsc_incr/*.mk')\
        +glob('*.mk')\
        +glob('umat/src/evpsc_incr/*.f')\
        +glob('umat/src/evpsc_incr/src/*.f')\
        +glob('umat/src/evpsc_incr/include/*.f')\
        +glob('src/vpsc_src/*.f')\
        +glob('src/vpsc_src/*.sub')\
        +glob('src/vpsc_src/*.dim')\
        +glob('pyscripts/*.py')

    fns.append('Makefile')
    fns.append('umat/src/evpsc_incr/Makefile')
    fns.append('pyscripts/archive.py')
    fns.append('pyscripts/statev_interface.py')
    fns.append('umat/src/evpsc_incr/diag.py')
    fns.append('umat/src/evpsc_incr/prebuild_vumat.py')
    fns.append('umat/src/evpsc_incr/lookup.py') # to run SA after FE
    fns.append('umat/src/evpsc_incr/abq_rst.py') # to run SA after FE
    fns.append('umat/src/evpsc_incr/gen_hist.py') # to run SA after FE
    fns.append('umat/src/evpsc_incr/plot_abq.py') # to run SA after FE

    tarfn=f'current_version_{hashcode}.tar'
    with tarfile.open(name=tarfn,mode='w') as myTarFile:
        for name in fns:
            myTarFile.add(name)
    ## move created tar
    shutil.move(tarfn,dst_base)
    ## extract the tar
    os.chdir(dst_base)
    with tarfile.open(name=os.path.join(dst_base,tarfn),mode='r') as myTarFile:
        myTarFile.extractall(filter='data')
    ## homing
    os.chdir(home)

def run_abq(dst,cpus,job='tmp',ilog=False,isolver=0,
            interactive=True,usermat=None):

    f"""
    Run abaqus using {ABQ}

    Arguments
    ---------
    dst
    cpus
    job
    ilog=False
    isolver=0
    interactive=True
    usermat=None ; name of the source code with the extension
    """
    import subprocess
    import sys
    from subprocess import PIPE
    home=os.getcwd()

    if type(usermat)==type(None):
        if isolver==0:
            usermat='umat.f'
        if isolver==1:
            usermat='vumat.f'

    ## test if the full pathed file is given (it should not be the case)
    upper,__src_file_name__=os.path.split(usermat)
    if len(upper)!=0:
        print('A single file name (without path) should be given ')
        print('for <usermat> in <run.run_abq> function')
        sys.exit(-1)

    ## for the case of umat (implicit), stripe of the file extension
    #if isolver==0:
    usermat=usermat.split('.')[0]

    try:
        os.chdir(dst)
        if isolver==0:
            if sys.platform=='win32':
                user=f'{usermat}-std.obj'
            elif sys.platform=='linux' or sys.platform=='linux2':
                user=f'{usermat}-std.o'
        elif isolver==1:
            user=f'{usermat}.f'

        cmd=[ABQ,f'job={job}',f'user={user}',
             'ask_delete=OFF','double',
             'mp_mode=threads','int',f'cpus={cpus}']
        if interactive:cmd.append('interactive')

        cmds=''
        for c in cmd: cmds=f'{cmds} {c}'
        print(f'cmd: {cmds}')

        if ilog:
            print('** with <run_stdout.log>')
            with open('run_stdout.log','w') as fo:
                proc=subprocess.run(cmd,stdout=fo)
        else:
            proc=subprocess.run(cmd,shell=True)
    except:
        print('--------------------------------------')
        print('** Error while attemping to run Abaqus')
        print('--------------------------------------')
        os.chdir(home)
        sys.exit(-1)
    os.chdir(home)

def write_evpsc_fe_in(probs,pathwork):
    """
    Arguments
    ---------
    probs
    pathwork
    """
    import sys
    import subprocess
    stv_dir=os.path.join(pathwork,'stv')
    cmd='python -u write_evpsc_fe_in.py --noel'

    probs=np.array(probs,dtype='int')
    noels=np.unique(probs[:,0])
    npts=np.unique(probs[:,1])

    for i in range(len(noels)):
        cmd=f'{cmd} {noels[i]}'
    cmd=f'{cmd} --npt'
    for i in range(len(npts)):
        cmd=f'{cmd} {npts[i]}'

    cmd=f'{cmd} --stv {os.path.join(pathwork,"stvdir")}'

    proc=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdo,stde=proc.communicate()

    if os.path.isdir(pathwork):
        shutil.move('evpsc_fe.in',pathwork)
    else:
        msg=f"** Error: could not find the given path: {path}"
        sys.stderr(msg)
        sys.exit(-1)


def gen_pathwork(tmp):
    """
    Given the <tmp> folder, generate a abaqus run-path
    in which the necessary source codes are to be located.
    We use <git rev-parse --short HEAD> to generate the
    hash code that is appended to the run-path.

    Argument
    --------
    tmp
    """
    import subprocess
    import os
    COD=datetime.datetime.now().strftime('%Y%m%d')
    CTT=datetime.datetime.now().strftime('%H%M%S')
    if os.name=='nt':
        # p=os.popen('git rev-parse --short HEAD')
        # hashcode=p.read().split('\n')[0]
        p=subprocess.Popen('git rev-parse --short HEAD'.split(),
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdo,stde=p.communicate()
        hashcode=stdo.decode().split('\n')[0]
        print(f'** hashcode obtained from <run.gen_pathwork>: {hashcode}')
    else:
        hashcode='None'
    pathwork=os.path.join(tmp,COD,CTT,hashcode,'umat','src','evpsc_incr')
    print('pathwork:',pathwork)
    os.makedirs(pathwork)
    return pathwork,hashcode


def add_and_select_material_card(inp,material_card):
    """
    Arguments
    ---------
    inp:           original abaqus input file,
           to which material_card is going to be added
    material_card: material card information for abaqus
    """
    import subprocess
    with open(inp,'r') as fo:
        lines=fo.read().split('\n')

    crit='*Material,'
    locs=[]
    for i in range(len(lines)):
        if lines[i][:len(crit)].lower()==crit.lower():
            locs.append(i)

    ## print("existing material cards:")
    mat_names=[]
    for i in range(len(locs)):
        loc=locs[i]
        mat_names.append(lines[loc].split('name=')[-1])
        print(f'Mat: {mat_names[i]}')

    ## search the line in which the material is assigned.
    line_ind_where_material_is_assigned=[]
    for i in range(len(lines)):
        for j in range(len(mat_names)):
            if f'material={mat_names[j].lower()}' in lines[i].replace(' ','').lower():
                line_ind_where_material_is_assigned.append(i)
                print(f'** Material is assigned: {lines[i]}')

    if len(line_ind_where_material_is_assigned)!=1:
        print('** Unexpected case: multiple lines exist where material is assigned')
        raise IOError('** Unexpected case')
        sys.exit(-1)

    ind=line_ind_where_material_is_assigned[0]
    print(f'before line: {lines[ind]}')
    head,tail=lines[ind].replace(' ','').split('material=')
    lines[ind]=f'{head}material=devpsc'
    print(f'after line: {lines[ind]}')

    #with open(inp,'r') as fo:
    #cnt_orig_inp=fo.read()
    #    lines=cnt_orig_inp.split('\n')
    cnt_orig_inp=''
    for l in lines:
        cnt_orig_inp=f'{cnt_orig_inp}{l}\n'

    ## insert new material card based on virutal dEVPSC simulation
    upper,lower=cnt_orig_inp.split(lines[locs[0]])
    with open(inp,'w') as fo:
        fo.write(f'{upper}{material_card}\n{lines[locs[0]]}{lower}')

    # sys.exit(0)

def pre_run_devpsc(iregime,conds,dt,edot,eps,pathwork,strfn=''):
    """
    If strfn (i.e., <STR_STR.OUT> file generated
    from dEVPSC-SA) is given, use it to generate <material data card>.

    Otherwise, conduct uniaxial (or bulge - not available yet)
    simulation with VPSC-SA or dEVPSC-SA and use the flow stress vs.
    strain curve directly from this result to generate
    the <material data card>.

    Arguments
    ---------
    iregime (0: VPSC, 1: delta EVPSC)
    ioptfs:   0: IJP 2020, 1: Still affine but with the secant
              compliance tensor; Applicable only for EVPSC
    conds
    dt
    edot
    eps
    pathwork
    strfn

    Returns
    -------
    mat_card  (material card required in Abaqus input file; <string> type)
    """

    import subprocess
    from MP.mat.mech import FlowCurve as FC
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import os

    if iregime==0:
        sys.stderr.write('*** Cannot use <pre_run_devpsc> with iregime=0 yet')
        sys.exit(-1)

    iodf=conds['iodf']

    print('*** Virtual experiment using dEVPSC to characterize a phenomenological model **')
    if conds['iodf']!=1:
        msg='**** cannot prepare a virtual characterization using dEVPSC unless iodf==1'
        sys.stderr.write(mgs)
        sys.exit(-1)
        raise IOError(msg)

    base=os.getcwd()
    if len(strfn)>0:
        # 0. In case <STR_STR.OUT> is given, use it to generate the material card
        myfc=FC()
        myfc.get_model(strfn)
    elif strfn=='':
        #--------------------------------------------------------------------------
        ## 1. Prep <base directories>
        print('**** prep base directory')
        for i in range(3): base=os.path.split(base)[0]
        dst_base=pathwork[::]
        for i in range(3): dst_base=os.path.split(dst_base)[0]
        #--------------------------------------------------------------------------
        ## 2. Build stand-alone code
        #..........................................................................
        # 2-1. parsers
        print('**** build parser')
        cmd=f'make parsers -j'
        p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=dst_base)
        stdo,stde=p.communicate()
        #..........................................................................
        # 2-2. Change ngrphmx and iodf
        print('**** change ngr, nphmx, iodf, and ioptfs')
        print(f"       nph: {conds['nph']}")

        cmd=f"python -u pyscripts/gen_dim.py --iregime {iregime} "
        cmd=f"{cmd} --ngrmx {conds['ntotgr']} --nphmx {conds['nph']} "
        cmd=f"{cmd} --fnsx "
        print('conds[sx]')
        print(conds['sx'])
        for i, _fn_ in enumerate(conds['sx']):
            cmd=f"{cmd} {os.path.join(os.getcwd(),_fn_)}"
        cmd=f"{cmd} -o {os.path.join(dst_base,'src/vpsc_src/parameters.dim')}"
        print(f'cmd: {cmd}')
        p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=dst_base)
        stdo,stde=p.communicate()
        print(stdo.decode())
        print(stde.decode())

        # iodf
        chg_parameter_val(
            fn=os.path.join(dst_base,'src/vpsc_src/evpsc_main.f'),
            parameter='iodf',val=conds['iodf'],iopt=0)
        # iodf
        chg_parameter_val(
            fn=os.path.join(dst_base,'src/vpsc_src/evpsc_main.f'),
            parameter='ioptfs',val=ioptfs,iopt=0)
        #..........................................................................
        # 2-3. Build using <Makefile> and <gfortran>
        cmd='make -j evpsc'
        print(f'**** <dst_base> in pre_run_devpsc: {dst_base} ')
        p=subprocess.Popen(cmd.split(),cwd=dst_base,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=p.communicate()
        try:
            print(stdout.decode(errors='ignore'))
            print(stderr.decode(errors='ignore'))
        except:
            print('** could not print out streams generated from running <make -j evpsc>')
            decoder1=chardet.detect(stdout)
            decoder2=chardet.detect(stderr)
            print(f'decoder1: {decoder1}')
            print(f'decoder2: {decoder2}')
            print(f'location: {dst_base}')

        #..........................................................................
        # 2-4. Restore the total number of grains for later FE simulations

        print(f'**** change dimension file with ntotgr')
        print(f"       nph: {conds['nph']}")

        #--------------------------------------------------------------------------
        ## 3. Generate input file for stand-alone run
        #..........................................................................
        # 3.1 generate the history - uniaxial tension along axis 1
        print(f'**** generate history file')
        histin=tempfile.mktemp(prefix="youngung_tmp_")
        cmd=f"python -u {os.path.join(base,'pyscripts','u_hist.py')} --ictrl 1 --edot {edot}"
        cmd=f"{cmd} --dt {dt} --eps {eps} --iopt 0 --fnout {histin}"
        p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        # with open(histin,'w') as fo: p=subprocess.Popen(cmd.split(),stdout=fo)
        #..........................................................................
        # 3.2 generate the process file
        print(f'**** generate process file')
        proc_cnt=f'** nproc\n1\n** ivgvar\n0\n{histin}'
        print(f'     - process content:\n{proc_cnt}')
        #..........................................................................
        # 3.3 use <vpsc7.in> for FE, and generate new one for stand-alone run.
        print(f'**** generate <vpsc7.in> file')
        shutil.copy(os.path.join(pathwork,'vpsc7.in'),dst_base)
        with open(os.path.join(dst_base,'vpsc7.in'),'r') as fo: cnt=fo.read()
        cnt=f'{cnt}{proc_cnt}'
        with open(os.path.join(dst_base,'vpsc7.in'),'w') as fo: fo.write(cnt)
        print(f'      - vpsc7.in file content: \n{cnt}')
        #--------------------------------------------------------------------------
        ## 4. Execute the stand-alone simulation
        print(f'**** Execution of stand-alone simulation')
        os.chdir(dst_base)
        if   sys.platform=='win32':                           fnexe='evpsc.exe'
        elif sys.platform=='darwin':                          fnexe='./evpsc'
        elif sys.platform=='linux' or sys.platform=='linux2': fnexe='./evpsc'
        p=subprocess.Popen([fnexe],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=dst_base)
        stdout,stderr=p.communicate()
        try:
            print(stdout.decode(errors='ignore'))
            print(stderr.decode(errors='ignore'))
        except:
            print('** could not print out streams generated from running EVPSC-SA')
            pass

        myfc=FC(); myfc.get_model('STR_STR.OUT')

    ## calculate elastic modulus
    fn_elsc=os.path.join(dst_base,'ELSC.OUT')
    elsc_row=np.loadtxt(fn_elsc)[0] ## Use only the initial ELSC
    #  0:1   - 1/aux3333
    #  1:37  - selav(i,j)
    # 37:41  - four bux3333
    # 41:41+81 - 81 compliance components
    # 41+81:41+81+81 - 81 stiffness components
    elsc_row=elsc_row[41+81:41+81+81]
    elsc=np.zeros((3,3,3,3))
    m=0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    elsc[i,j,k,l]=elsc_row[m]
                    m=m+1
    fn=os.path.join(dst_base,'aniso_elsc.out')
    with open(fn,'w') as fo:
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        fo.write('%12s'%f'E_{i+1}{j+1}{k+1}{l+1}')
        fo.write('\n')
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        fo.write('%12.1f'%elsc[i,j,k,l])


    ijkl='1111 1122 2222 1133 2233 3333 1112 2212 3312 1212 1113 2213 3313 1213 1313 1123 2223 3323 1223 1323 2323'
    ijkl=ijkl.split()
    indices=[]
    for m, ijkl in enumerate(ijkl):
        i,j,k,l=ijkl
        i,j,k,l=map(int,[i,j,k,l])
        indices.append([i-1,j-1,k-1,l-1])
    elast_card='*Elastic, type=ANISOTROPIC\n'
    for m, ijkl in enumerate(indices):
        i,j,k,l=ijkl
        val=elsc[i,j,k,l]
        elast_card=f'{elast_card} %9.0f,'%val

    ## this algorithm might need refinement
    ys0=250.
    dys=ys0/25. #.10

    x,y=np.array(myfc.epsilon_vm),np.array(myfc.sigma_vm)
    fig=plt.figure();ax=fig.add_subplot()
    ax.plot(x,y,label='Total strain vs. stress')

    # resi=1e-10
    # k=0
    # while resi<1e-4:
    flt=(y>0)*(y<ys0)
    X=x[flt]
    Y=y[flt]
    rst=np.polyfit(X,Y,deg=1,full=True)
    Z,resi=rst[0],rst[1]
    ys0=ys0+dys
    # k=k+1
    # if (k>1000):
    #     print("something wrong")
    #     sys.exit(-1)

    ## determine elastic modulus, poisson's ratio, and yield stress
    ys0=ys0-dys
    flt=(y>0)*(y<ys0)

    e33=myfc.epsilon[2,2,flt]
    e22=myfc.epsilon[1,1,flt]
    e11=myfc.epsilon[0,0,flt]
    nu1=-e33[-1]/e11[-1]
    nu2=-e22[-1]/e11[-1]
    nu=(nu1+nu2)*0.5 ## use an average of nu1 and nu2

    if nu>0.495:
        sys.stderr.write(f"** Abaqus does not like poisson ratio exceeding 0.495")
        sys.stderr.write(f"** The calculated nu is {nu} but is reset to 0.49")
        nu=0.49


    emod=(Y[-1]-Y[0])/(X[-1]-X[0])

    epl=x-y/emod
    flt=epl>0.002 # offset strain
    epl=epl[flt]
    sig=y[flt]


    flg=np.diff(sig)>0
    sig=sig[1:][flg]
    epl=epl[1:][flg]


    ax.plot(epl,sig,label='Plastic strain vs stress (not shifted)')
    epl=epl-epl[0]
    ax.plot(epl,sig,label='Plastic strain vs stress (zeroed)')
    ax.legend()

    fig.savefig('strstr_devpsc.pdf')
    os.chdir(base)
    dat=np.array([x,y])

    ## write down material card
    mat_card=f'*Material, name=devpsc\n*Elastic\n{emod}, {nu}'
    mat_card=f'{mat_card}\n*Plastic'
    ## remove the last one (as it pertains to repeat but without update in sa run)
    for i in range(len(epl)-1):
        # mat_card=f'{mat_card}\n{sig[i]}, {epl[i]}'
        mat_card=f'{mat_card}\n %15.7e, %15.7e'%(sig[i],epl[i])

    return mat_card, elast_card

def check_prerequisites():
    """
    check if intel fortran compiler <ifort> is available
    and is running in win32
    """
    import sys
    import subprocess
    cmd=['ifort']
    kwargs=dict(stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p=subprocess.Popen(cmd,**kwargs)
    stdo,stde=p.communicate()
    lines=stde.decode(errors='ignore').split('\n')

    willIrun=True
    if sys.platform=='win32':
        if not('Intel' in lines[1]):
            willIrun=False
    else:
        sys.stderr.write('** Warning: The check prerequisite test is not ready for system other than windows\n')

    return willIrun

if __name__=='__main__':
    from pathlib import Path
    import argparse
    import sys
    import tempfile
    import subprocess
    import os
    from build_umat import build_prep, build_sa, build_fe
    dir_init=os.getcwd()
    # --------------------------------------------------------------------------------------------
    ## check if the prerequisites are available
    willIrun=check_prerequisites()
    if not(willIrun):
        sys.stderr.write('**Error, I do not expect the system will run - check ifortran/abaqus/environment\n')
        sys.exit(-1)
    # --------------------------------------------------------------------------------------------
    sys.stdout.write("""--------------------------------------------------------------------------\n""")
    sys.stdout.write("""Run Abaqus/standard or explicit simulation using dEVPSC-FE umat (or vumat)\n""")
    sys.stdout.write('\n')
    sys.stdout.write("""  This script allows to neatly run dEVPSC-FE with prescribed conditions\n""")
    sys.stdout.write("""  saved in the given files that include necessary conditions. \n""")
    sys.stdout.write("""    - Explicit version is under progress...\n""")
    sys.stdout.write("""          Youngung Jeong; yjeong@changwon.ac.kr\n""")
    sys.stdout.write("""--------------------------------------------------------------------------\n\n""")
    # --------------------------------------------------------------------------------------------
    ## 1. Argument parsing
    ### 1.1 List of arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('--cpus',type=int,help='Total number of CPUs',default=1)

    parser.add_argument('--fn',type=str,help='Input file to <run.py>; examples are relocated in <umat/src/evpsc_incr/tests> folder.',default='umat.in')
    parser.add_argument('--iregime',type=int,default=2,help='1: VPSC; 2: delta-EVPSC (incremental formulation)')
    parser.add_argument('--ioptfs',type=int,default=-1,help='0: IJP 2020, 1: Still affine but with a secant compliance tensor; Applicable only for EVPSC')

    parser.add_argument('--iodf',type=int,help='0: discrete; 1: odf (mtex), 2: random, 3: EBSD, 4: Euichan, 5: shuffle')

    parser.add_argument('--ngrs',nargs='+',default=[0],help='The total number of grains per each phase within each FE integration point.')

    parser.add_argument('--idd', action='store_true',help='Boolean flag to switch on the DD hardening model; switch in <statev_io.f>')
    parser.add_argument('--irgvb', action='store_true',help='Boolean flag to switch on the RGVB-DD hardening model; switch in <statev_io.f>')

    parser.add_argument('--pathwork',type=str,help='Name of path, in which FE is to run. If not given, it is generated within <run.py>.')

    parser.add_argument('--j2',action='store_true',help='Run with J2 theory with utilizing the flow curve generated from stand-alone EVPSC simulation')
    parser.add_argument('--j2strfn',type=str,help='STR_STR.OUT file if exists (optional); instead of actually running stand-alone EVPSC simulation, use the given <STR_STR.OUT> file to determine the flow stress vs. strain curve',default='')
    parser.add_argument('--isa', action='store_true',help='Build and run SA after FE simulation')
    parser.add_argument('--instance',type=str,help='Abaqus instance name to check from <lookup.py> via <abq_rst.py>',default='part-1-1')
    parser.add_argument('--inp',type=str,help='Abaqus input file')#,default='')
    parser.add_argument('--elastmod',action='store_true',help='Calculate elastic modulus (with 10,000 grains per each phase) and exit')
    parser.add_argument('--usermat',type=str,default='',help='Your user material subroutine; which is set to be <umat.f> unless is specified.')

    parser.add_argument('--EBSDrecipe',type=int,help='1: C-ring; 2-6: POSCO welds')
    parser.add_argument('--fntx',nargs='+',help='Replacing texture(s)')

    parser.add_argument('--idebug', action='store_true',help='Boolean flag to switch on debugging options')

    args=parser.parse_args()
    args.ngrs=np.array(args.ngrs,dtype='int')

    ### 1.2 Testing arguments
    if args.iregime==2 and args.ioptfs==-1:
        sys.stderr.write('** For EVPSC, (iregime==1), ioptfs should be specified')
        sys.exit(-1)

    # --------------------------------------------------------------------------------------------
    ## 2. Read base conditions from the given file 'args.fn',
    ## which can be overriden by arguments parsed.
    sys.stdout.write( "** Read input file to <run.py>\n")
    sys.stdout.write(f'*** args.fn :{args.fn}\n')
    conds=read_in(fn=args.fn)
    #    * 2-1. Specify the number of grains
    ntotgr=0 ## number of grains pertaining to each integration point
    # if not for elastic modulus calculations
    if not(args.elastmod):
        for iph in range(conds['nph']): ntotgr=ntotgr+args.ngrs[iph]
    if args.elastmod:      # if for elastic modulus calculations, use 10000 grains per each phase
        for iph in range(conds['nph']): ntotgr=ntotgr+10000
    if conds['ipx']!=0  and 0 in args.ngrs:
        sys.stderr.write("** Warning: make sure you are using the same ngrs for px and FE runs\n")
    if conds['iodf']!=0 and 0 in args.ngrs:
        sys.stderr.write('** Error: When odf/EBSD is given (one must provide <ngrs>\n')
        sys.exit(-1)

    conds.update(ntotgr=ntotgr)

    ## - replace the texture
    if type(args.fntx)!=type(None):
        conds.update(text=args.fntx)
    if type(args.iodf)!=type(None):
        conds.update(iodf=args.iodf)

    # --------------------------------------------------------------------------------------------
    ## 3. Create temporary folder (force to prefer d:\\temp as a default.)
    sys.stdout.write("** Create temporary folder\n")
    if sys.platform=='win32':
        ## for Windows server 2019 at the lab, the one next to the DIC system.
        tmp='d:\\temp'
        if not(os.path.isdir(tmp)): tmp='c:\\temp'
    elif sys.platform=='darwin': ## Mac OS
        tmp='/tmp'
    else:
        user=os.environ['USER']
        tmp=os.path.join('/scratch',user)
        sys.stderr.write(f'** Unexpected system platform {sys.platform}\n')
        sys.stderr.write(f'tmp: {tmp}\n')
        sys.stderr.write('** exit\n')
        sys.exit(-1)

    # --------------------------------------------------------------------------------------------
    ## 4. Create <pathwork> (and 'hashcode' of the current commit)
    #     or use the given one by dumping auto-generated one
    sys.stdout.write("** Create <pathwork>\n")
    if type(args.pathwork)==type(None):
        pathwork, hashcode=gen_pathwork(tmp)
    else: # in case pathwork is not given
        if not(os.path.isdir(args.pathwork)):
            sys.stderr.write(f'** Could not find pathwork {args.pathwork}\n')
        else:
            sys.stdout.write('** Given path work exists, remove and recreate it\n')
            shutil.rmtree(args.pathwork)
        pathwork=os.path.join(args.pathwork,'umat','src','evpsc_incr')
        Path(pathwork).mkdir(parents=True,exist_ok=True)
        # <_pathwork_> from gen_pathwork is not used.
        # Only the hashcode generated from gen_pathwork is used
        # for naming the tar file to be created from
        _pathwork_,hashcode=gen_pathwork(tmp)
        ## copy any dim file

    with open('run.log',mode='a') as fo: fo.write(f'\n{pathwork}')
    shutil.copy(args.fn, pathwork)

    dst_base=pathwork[::]
    for i in range(3): dst_base=os.path.split(dst_base)[0]

    # --------------------------------------------------------------------------------------------
    ## 5. Copy the necessary abaqus environment file when Nurion cluster is used.
    sys.stdout.write("** Copy necessary envir file\n")
    if   sys.platform=='win32':  pass
    elif sys.platform=='darwin': pass
    elif sys.platform=='linux' or sys.platform=='linux2':
        shutil.copy('/apps/commercial/abaqus/2020/EstProducts/2020/linux_a64/SMA/site/abaqus_v6.env',pathwork)

    # --------------------------------------------------------------------------------------------
    ## 6. Prep <vpsc7.in> in which VPSC numerical conditions are given
    #   by using "gen_inp.py" module located three directory-levels below.
    sys.stdout.write("** Prep vpsc7.in file\n")
    ref=os.path.join(pathwork,'vpsc_ref.inp')
    cmd=f"python -u {os.path.join('pyscripts','gen_inp.py')} --fnout {ref} --nph {conds['nph']} --wph"
    for iph in range(conds['nph']): cmd=f"{cmd} {conds['wph'][iph]}"
    sys.stdout.write(f'cmd: {cmd}\n')
    p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd='../../../')
    stdo,stde=p.communicate()
    #  Create the vpsc7.in file by 'replacing' place-holding lines by the wanted input lines.
    sys.stdout.write(f"conds['text']: {conds['text']}\n")
    prep_vpscin_nph(dst=pathwork,idiffs=conds['idiff'],textures=conds['text'],
                    sxs=conds['sx'],diffs=conds['diff'],
                    ngrs=args.ngrs,iodf=conds['iodf'],ref=ref)

    # --------------------------------------------------------------------------------------------
    ## 7. Use write <evpsc_fe.in>
    sys.stdout.write("** Write <evpsc_fe.in>\n")
    if len(conds['probs'])>0:
        write_evpsc_fe_in(conds['probs'],pathwork)

    # --------------------------------------------------------------------------------------------
    ## 8. Archive
    # Create and extract tar in the targeted destination, i.e., <pathwork>
    # It changes <dim> files belonging to <pathwork>
    sys.stdout.write("** Archive\n")
    archive(hashcode,pathwork,ntotgr)
    #     - Copy Abaqus input file to destination as <tmp.inp>
    job_name='tmp' ## the default job name is 'tmp', which can be changed later.
    if type(args.inp)==type(None):
        fnabqinput_origin=conds['abqinp']
    else:
        fnabqinput_origin=args.inp

    check_abqinp(fnabqinput_origin,conds['isolver'])
    shutil.copy(fnabqinput_origin,os.path.join(pathwork,f'{job_name}.inp'))
    #      - Copy the <args.usermat> file  to the pathwork, if given.
    sys.stdout.write(f"conds['isolver'] : {conds['isolver']}\n")
    if args.usermat!='':
        if conds['isolver']==0: dst_fn=os.path.join(pathwork,'umat.f')
        elif conds['isolver']==1: dst_fn=os.path.join(pathwork,'vumat.f')
        else:
            sys.stderr.write('** Error: unexpected isolver\n')
            sys.exit(-1)
        shutil.copy(args.usermat,dst_fn)

    # --------------------------------------------------------------------------------------------
    ## 9-1. Change ipxx/iodf flags in <sdvini.f>;
    ## This subroutine is required only for UMAT not VUMAT
    sys.stdout.write("** Change ipxx/iodf in <sdvini.f>, <uexternaldb.f>\n")
    if conds['isolver']==0:
        # chg_ipx(conds['ipx'],  os.path.join(pathwork,'sdvini.f')) # this might be revised to 'ipxx' not 'ipx'
        chg_parameter_val(fn=os.path.join(pathwork,'sdvini.f'),parameter='ipxx',val=conds['ipx'],iopt=0)
        # chg_iodf(conds['iodf'],os.path.join(pathwork,'sdvini.f'))
        chg_parameter_val(fn=os.path.join(pathwork,'sdvini.f'),parameter='iodf',val=conds['iodf'],iopt=0)
        # chg_iodf(conds['iodf'],os.path.join(pathwork,'src','uexternaldb.f'))
        chg_parameter_val(fn=os.path.join(pathwork,'src','uexternaldb.f'),parameter='iodf',val=conds['iodf'],iopt=0)
        #chg_ebsdfn(os.path.join(pathwork,'src','uexternaldb.f'),ref='parameter(ebsd_filename=',ebsd_fn='texture_0')
        chg_parameter_val(fn=os.path.join(pathwork,'src','uexternaldb.f'),parameter='ebsd_filename',val='texture_0',iopt=1)

        sys.stdout.write(f"** Check: {os.path.join(pathwork,'src','uexternaldb.f')}\n")

        if conds['iodf']==3: ## applies only when using *.ang EBSD data file.
            sys.stdout.write("** Change nrowmx in <data_grain_ebsd.f>, <evpsc_input_fe.f>, \n")
            sys.stdout.write("  <select_points_from_ebsd.f>, <sdvini.f>, <uexternaldb.f> and <read_ebsd.f>\n")
            with open(os.path.join(pathwork,'texture_0'),'r') as fo:
                nrowmx=len(fo.readlines()) ## this contains the possible number of head (comment) lines
            files_including_nrowmx=[os.path.join(dst_base,'src','vpsc_src','data_grain_ebsd.f'),
                                    os.path.join(dst_base,'src','vpsc_src','evpsc_input_fe.f'),
                                    os.path.join(dst_base,'src','vpsc_src','select_points_from_ebsd.f'),
                                    os.path.join(dst_base,'src','vpsc_src','read_ebsd.f'),
                                    os.path.join(pathwork,'sdvini.f'),
                                    os.path.join(pathwork,'src','uexternaldb.f')]
            for i, fn in enumerate(files_including_nrowmx):
                chg_parameter_val(fn,'nrowmx',nrowmx)

    # --------------------------------------------------------------------------------------------
    ## 9-2. Change EBSD recipe in <select_points_from_ebsd.f>
    if conds['iodf']==3 and type(args.EBSDrecipe)==type(None):
        sys.stderr.write("** Error: iodf is 3 but no specific ")
        sys.stderr.write(" <EBSDrecipe> is given to <run.py>\n")
        sys.exit(-1)
    if type(args.EBSDrecipe)!=type(None):
        fnsrc='select_points_from_ebsd.f'
        sys.stdout.write(f"** Select EBSD recipe option in <{fnsrc}>\n")
        fnsrc=os.path.join(dst_base,'src','vpsc_src',fnsrc)
        chg_parameter_val(fnsrc,'irecipe',args.EBSDrecipe)

    # --------------------------------------------------------------------------------------------
    ## 10. If plans to build & run Stand-alone code, do so.
    sys.stdout.write("** Build & run SA simulation if necessary\n")
    if args.isa or args.elastmod:
        sys.stdout.write("*** change iodf\n")
        # chg_iodf(conds['iodf'],os.path.join(pathwork,'../../../src/vpsc_src/evpsc_main.f'))
        chg_parameter_val(fn=os.path.join(pathwork,'../../../src/vpsc_src/evpsc_main.f'),parameter='iodf',val=conds['iodf'],iopt=0)

    # --------------------------------------------------------------------------------------------
    ## 11-1. If plans to run J2 plasticity theory run (based on von Mises)
    if args.j2:
        sys.stdout.write("*** J2 option activated, run SA simulation\n")
        #   11-1. If j2strfn is given, use it to generate the material card,
        #   11-2. Else, run dEVPSC to generate the material card based on
        #      "virtual" test by simulating uniaxial tension. Use the resulting
        #      stress-strain curve data file "STR_STR.OUT"
        mat_card, elast_card =pre_run_devpsc(
            iregime=args.iregime,ioptfs=args.ioptfs,conds=conds,dt=0.5,edot=0.001,eps=0.20,
            pathwork=pathwork,strfn=args.j2strfn)
        ## create new job file name <virtual_devpsc.inp>
        job_name='virtual_devpsc' # change the job name!
        jobfn_devpsc=os.path.join(pathwork,f'{job_name}.inp')
        shutil.copy(os.path.join(pathwork,'tmp.inp'), jobfn_devpsc)
        ## add the above-generated material card to newly created input file
        sys.stdout.write("*** Putting J2 result to material card")
        sys.stdout.write("in the given ABQ inp file\n")
        add_and_select_material_card(inp=jobfn_devpsc,material_card=mat_card)

    ## 11-2, If elastic modulus is used.
    if args.elastmod:
        sys.stdout.write("*** if elastic card is requested\n")
        mat_card, elast_card = pre_run_devpsc(
            iregime=args.iregime,conds=conds,dt=0.1,edot=0.001,eps=0.001,
            pathwork=pathwork,strfn=args.j2strfn)
        sys.stdout.write('\n\n**********************************\n')
        sys.stdout.write('** elastic card obtained as below:\n')
        sys.stdout.write('**********************************\n\n')
        sys.stdout.write(elast_card)

    # --------------------------------------------------------------------------------------------
    ## 12. Build EVPSC-SA and umat-std.obj or attempt to build single precision vumat for checking
    sys.stdout.write("** Build EVPSC-SA & umat-std.obj\n")
    if conds['isolver']==0 or True:
        sys.stdout.write('** Build begins\n')
        sys.stdout.write(f"conds['iodf']: {conds['iodf']}\n")
        sys.stdout.write(f"conds['sx']: {conds['sx']}\n")
        sys.stdout.write(f"iregime: {args.iregime}\n")
        sys.stdout.write("*** Build prep\n")

        build_prep(dst=pathwork,conds=conds,
                   ioptfs=args.ioptfs,
                   idebug=args.idebug,
                   irgvb=args.irgvb,idd=args.idd,
                   iregime=args.iregime)
        if args.isa:
            sys.stdout.write("*** Build SA\n")
            build_sa(dst=dst_base)
        sys.stdout.write("*** Build FE\n")
        build_fe(dst=pathwork,isolver=conds['isolver'])

    # --------------------------------------------------------------------------------------------
    ## 13. if statev_px is present (maybe this is valid only for standard/implicit runs)
    if conds['ipx']==1 or conds['ipx']==2:
        sys.stdout.write("** Consder statev_px file\n")
        shutil.copy(conds['pxfn'],os.path.join(pathwork,'statev_px'))

    # --------------------------------------------------------------------------------------------
    ## 14. run Abaqus
    sys.stdout.write("** Abaqus simulation begins\n")
    run_abq(dst=pathwork,cpus=args.cpus,ilog=False,isolver=conds['isolver'],
            interactive=True,job=job_name)

    # --------------------------------------------------------------------------------------------
    ## 15. Check up in comparison with stand-alone case.
    # run <lookup.py> - running stand-alone code to re-simulate by using the history of a chosen
    # material point for checking up.
    if len(conds['probs'])>0 and args.isa and not(args.j2):
        sys.stdout.write("** Look up\n")

        ihardlaw=0 # voce
        if args.idd:
            ihardlaw=20 # DD (B&T)
        if args.irgvb:
            ihardlaw=23 # DD (RGVB)

        for i, pair in enumerate(conds['probs']):
            noel, npt = pair
            cmd=f"python -u lookup.py -p {pathwork}"
            cmd=f"{cmd} -noel {noel} -npt {npt}"
            cmd=f"{cmd} -instance_name {args.instance} -ihardlaw {ihardlaw}"
            sys.stdout.write(f'cmd: {cmd}\n')

            kwargs=dict(stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=dir_init)

            p=subprocess.Popen(cmd.split(),**kwargs)
            stdo,stde=p.communicate()
            print(stdo.decode())
            print(stde.decode())
