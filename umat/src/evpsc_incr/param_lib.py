import numpy as np
import matplotlib.pyplot as plt
import subprocess,shutil

import datetime
import os
import tempfile
import pickle
from datetime import datetime,date
from time import strptime


def calc_secs_explicit(h,m,s,day,year,ampm):
    # convert 12 hour clock to 24 hour clock
    # print('s in exp:',s)
    if ampm=='AM' and h==12:
        h=0
    elif ampm=='PM' and h==12:
        h=0
    if ampm=='PM':
        h=h+12

    secs=h*3600+m*60+s #+ampm
    # print('secs:',secs)
    # print('mins:',secs/60.)
    # print('hour:',secs/3600.)
    # print('-------------------')

    return secs

def calc_secs(l):
    """
    Year,Months, days are not accounted (required to use datatime class)
    """
    # print('** l in calc_secs:',l) ##  Mon 18 Jan 2021 10:32:47 AM KST
    ln=len(l.split())
    # print('ln:',ln)

    if ln==3:
        mdy,time,ampm=l.split()
        mon,day,year=mdy.split('/')
        h,m,s=time.split(':')
    elif ln==7:
        # print('l in calc_secs:',l)
        'Mon 18 Jan 2021 10:32:47 AM KST'
        date,day,mon,year,time,ampm,dum=l.split()
        h,m,s=time.split(':')
    else:
        print('something went wrong')
        print('l:',l)
        raise IOError

    #print('h,m,s:',h,m,s)
    h=int(h)
    m=int(m)

    # convert 12 hour clock to 24 hour clock

    if ampm=='AM' and h==12:
        h=0
    elif ampm=='PM' and h==12:
        h=0

    if ampm=='PM':
        h=h+12

    #print('ampm:',ampm)
    secs=h*3600+m*60+float(s) #+ampm

    # print('secs:',secs)
    # print('mins:',secs/60.)
    # print('hour:',secs/3600.)
    # print('-------------------')

    return secs
def time_delta(l1,l2):
    # print('l1:',l1)
    # print('l2:',l2)
    dsec=calc_secs(l2)-calc_secs(l1)
    # print('------------------')
    # print('dsec :',dsec)
    # print('dmin :',dsec/60.)
    # print('dhour:',dsec/3600.)
    # print('------------------')


    try: # windows
        mdy,time,ampm=l1.split()
        mon,day,year=mdy.split('/')
        year=int(year)
        mon=int(mon)
        day=int(day)
    except:
        date_,day,mon,year,time,ampm,dum=l1.split()
        year=int(year)
        day=int(day)
        mon=strptime(mon,'%b').tm_mon
        mon=int(mon)

    date1=date(year,mon,day)

    try: # windows
        mdy,time,ampm=l2.split()
        mon,day,year=mdy.split('/')
        year=int(year)
        mon=int(mon)
        day=int(day)
    except:
        date_,day,mon,year,time,ampm,dum=l2.split()
        year=int(year)
        day=int(day)
        mon=strptime(mon,'%b').tm_mon
        mon=int(mon)

    date2=date(year,mon,day)

    date_delta = date2-date1
    date_delta = date_delta.total_seconds()


    # print('------------------')
    # print('date_delta:',date_delta/3600.)
    # print('------------------')


    return dsec+date_delta
def gen_rand_tex_fn(ngr=100):
    from TX import cmb
    gr=cmb.random(phi1=360,phi2=360,phi=360,ngrain=ngr,iplot=False)
    fid,fn=tempfile.mkstemp()
    with open(fn,'w') as fo:
        fo.write('dum\ndum\ndum\n B %i\n'%ngr)
        for i in range(len(gr)):
            fo.write('%e %e %e %e\n'%(gr[i,0],gr[i,1],gr[i,2],gr[i,3]))
    fn =fn.replace('\\','/')
    print('rand_texture:',fn)
    return fn

def replace_rand_texture_in_input(ngr):
    evpsc_input="""1                          number of integration points per an element (should be hardwired)
1                          number of phases per each integration point
1.0                        relative vol. fract. of phases (wph(i))
*INFORMATION ABOUT PHASE #1
0   0   25                    grain shape contrl, fragmentn, crit aspect ratio
1.0  1.0  1.0                 initial ellipsoid ratios (dummy if ishape=4)
0.0  0.0  0.0                 init Eul ang ellips axes (dummy if ishape=3,4)
* name and path of texture file (filetext);
%s
* name and path of single crystal file (filecrys)
c:/users/user/repo/evpsc-yj-local/umat/src/evpsc_incr/matData/lc/lc.txt
* name and path of grain shape file (dummy if ishape=0) (fileaxes)
dummy
* diffraction (0: none,  1: LANL   2: NIST)
0
c:/users/user/repo/evpsc-yj-local/examples/bench_MMM2018/fcc.dif
*PRECISION SETTINGS FOR CONVERGENCE PROCEDURES (default values)
0.001 0.001 0.001 0.01  ## ERRS,ERRD,ERRM,ERRSO
50  50  15    itmax:  max # of iter, external, internal and SO loops
0  2 10 2   irsvar & xrsini,xrsfin,xrstep (dummy if irsvar=0)
1              ibcinv (0: don't use <Bc>**-1, 1: use <Bc>**-1 in SC eq)
*INPUT/OUTPUT SETTINGS FOR THE RUN (default is zero)
0              irecover:read grain states from POSTMORT.IN (1) or not (0)?
0              isave:   write grain states in POSTMORT.OUT at step 'isave'?
0              icubcomp:calculate fcc rolling components?
0              nwrite (frequency of texture downloads)
*MODELING CONDITIONS FOR THE RUN
1              interaction (0:FC,1:affine,2:secant,3:neff=10,4:tangent,5:SO)
1  1  1        iupdate: update orient, grain shape, hardening
0              nneigh (0 for no neighbors, 1 for pairs, etc.)
0              iflu (0: don't calc, 1: calc fluctuations)
"""
    fntex=gen_rand_tex_fn(ngr)
    evpsc_input=evpsc_input%fntex

    with open('vpsc7.in','w') as fo:
        fo.write(evpsc_input)

def find_param_line(fn='../../../src/vpsc_src/parameters.dim',param='NGRMX'):
    """
    Find the line number correspoonding to the parameter statement for the
    given keyword.

    Arguments
    ---------
    fn
    param

    Returns
    -------
    i  : line number
    lines : the original files read in lines
    """
    with open(fn,'r') as fo:
        lines=fo.readlines()

    for i in range(len(lines)):
        print(lines[i])
    ## find line that contains NGRMX
    for i in range(len(lines)):
        try:
            item=lines[i].split('=')[0]
        except:pass
        else:
            if '=' in lines[i]:
                if param in lines[i].upper():
                    return i,lines
    raise IOError(f'Could not find {kwd} line')

def replace_param(fn,param,ngr):
    """
    Replace parameter statement

    Arguments
    ---------
    fn
    param
    ngr
    """
    ind,lines=find_param_line(fn,param)

    if int(lines[ind].split('=')[-1].split(')')[0])==ngr:
        print('** No need to replace')
        return

    lines[ind]=f'      PARAMETER({param}={ngr})\n'

    with open(fn,'w') as fo:
        for i in range(len(lines)):
            fo.write('%s'%lines[i] )

def compile_umat_param(fn,param,ival):
    """
    1. Replace the parameter <param> in file <parameters.dim> with new value <ival>
    2. compile to obtain umat-std.obj and standardU.dll
    3. Copy the obtained object and dll files to different files and return the file names


    Arguments
    ---------
    fn         ; EVPSC-UMAT's incr_evpsc_fe.dim file
    param      ; Parameter to be modified in <incr_evpsc_fe.dim>
    ival       ; New parameter for the given <param> keyword.

    Returns
    -------
    newumatfn  ; new umat file name
    newdllfn   ; new dll file name
    """
    import time
    replace_param(fn,param,ival)
    cmd=['make','-i']
    t0=time.time()
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    t1=time.time()

    if t1-t0 < 10 or True:
        print('--------------------------------------------------')
        print('Elapsed time for compile/Linking:',t1-t0,' [secs]')
        print('--------------------------------------------------')
    # else:
    #     stdo,stde=p.communicate()
    #     print('stdo')
    #     print(stdo)
    #     stdo=stdo.decode('utf-8')
    #     lines=stdo.split('\n')

    #     i0=-1
    #     i1=0
    #     for i in range(len(lines)):
    #     # print(i,'--',lines[i])
    #         if 'Begin Compiling Abaqus' in lines[i]:
    #             i0=i+1
    #         if 'End Linking Abaqus' in lines[i]:
    #             i1=i+1

    #     dtime=time_delta(lines[i0],lines[i1])
    #     print('--------------------------------------------------')
    #     print('Elapsed time for compile/Linking:',dtime,' [secs]')
    #     print('--------------------------------------------------')

    ## save umat-std.obj and standardU.dll
    newumatfn='umat-std_%s_%i.obj'%(param,ival)
    newdllfn='standardU_%s_%i.dll'%(param,ival)
    shutil.move('umat-std.obj',newumatfn)
    shutil.move('standardU.dll',newdllfn)

    print('newumat:',newumatfn)
    print('newdllfn:',newdllfn)

    return newumatfn, newdllfn

def run_abq(inp,mem,cpu,ngrmx,umat,dll,test):
    """
    Arguments
    ---------
    inp: abaqus input file name (without extension 'inp')
    mem: memory (string, e.g., '5 gb')
    cpu: # of cpus to use (integer)
    ngrmx: Maximum # of allowed grains (integer)
    umat: umat objective file to use
    dll: standardU.dll file to use
    test: if True, not running Abaqus.

    Returns
    -------
    time
    """
    ## copy obj and dll files.
    print('umat: %s \ndll: %s'%(umat,dll))
    if umat!='umat-std.obj': shutil.copy(umat,'umat-std.obj')
    if dll!='standardU.dll': shutil.copy(dll,'standardU.dll')

    ## make sure we have random texturew with a correct number of grains
    replace_rand_texture_in_input(ngrmx)

    ## Ensure inp has a correct DEPVAR
    replace_param('../../../src/vpsc_src/parameters.dim','NGRMX',ngrmx)

    cmd=['make','-i','-j','test']
    p=subprocess.Popen(cmd)#,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    cmd=['test']
    p=subprocess.Popen(cmd)#,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()

    cmd=['C:/SIMULIA/EstProducts/2020/win_b64/code/bin/ABQLauncher.exe',
         'interactive',
         'job=%s'%inp,
         'user=umat-std.obj',
         'ask_delete=OFF',
         'cpus=%i'%cpu,
         'memory="%s"'%mem]

    cmds=cmd[0]
    for s in cmd[1:]:
        cmds='%s %s'%(cmds,s)

    print(cmds)

    if test: return np.nan

    p=os.popen(cmds)
    cnt=p.read()

    lines=cnt.split('\n')
    i0=-1
    i1=0
    for i in range(len(lines)):
        if 'Run standard.exe' in lines[i]:
            print('indx:',i,lines[i])
            print(lines[i-1])
            i0=i-1

        if 'End Abaqus/Standard Analysis' in lines[i]:
            print('indx:',i,lines[i])
            print(lines[i-1])
            i1=i-1

    print('i0:',i0,' i1:',i1)

    try:
        dtime=time_delta(lines[i0],lines[i1])#/3600.
        print('-------------------------------------------')
        for i in range(len(lines)):
            print(lines[i])
        print('-------------------------------------------')
    except:
        print('Something must have gone wrong')
        print('Below is the actual lines from standard.exe')
        print('-------------------------------------------')
        for i in range(len(lines)):
            print(lines[i])
        print('-------------------------------------------')
        return np.nan

    return dtime

# In[19]:


if __name__=='__main__':

    mems=['2gb','3gb']
    cpus=[4,6,8]
    cmd=['C:/SIMULIA/EstProducts/2020/win_b64/code/bin/ABQLauncher.exe',
         'interactive',
         'job=UMAT_DEV_VEL',
         'user=umat-std.obj',
         'ask_delete=OFF',
         'cpus=8']
    cmds=''
    for i in range(len(cmd)):
        cmds='%s %s'%(cmds,cmd[i])
        print(cmds)

    replace_param('../../../src/vpsc_src/parameters.dim','NGRMX',50)

    inps=['UMAT_DEV_VEL_1','UMAT_DEV_VEL_8','UMAT_DEV_VEL_27']#,'UMAT_DEV_VEL_64']

    delta_times=np.zeros((len(mems),len(cpus),len(ngrmxs),len(inps)))

    for i in range(len(mems)):
        for j in range(len(cpus)):
            for k in range(len(ngrmxs)):
                for l in range(len(inps)):
                    delta_times[i,j,k,l]=run_abq(inp=inps[l],mem=mems[i],cpu=cpus[j],
                                                 ngrmx=ngrmxs[k],umat=umats[k],dll=dlls[k],test=False)

    with open('param_results.pck','w') as fo:
        pickle.dump(fo,mems,cpus,ngrmxs,delta_times)
