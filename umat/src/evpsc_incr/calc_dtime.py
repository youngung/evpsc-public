import numpy as np
from param_lib import time_delta, calc_secs_explicit, calc_secs


def read_line_explicit_windows(line):
    month,day,year=map(int,line.split()[0].split('/'))
    h,m,s=map(int,line.split()[1].split(':'))
    ampm=line.split()[2][:2]
    # print(year,month,day,h,m,s)
    s0=calc_secs_explicit(h,m,s,day,year,ampm)
    return s0


def read_line_explicit_linux(line):
    # 'Mon 18 Jan 2021 10:32:47 AM KST'
    date,day,mon,year,time,ampm,dum=line.split()
    h,m,s=time.split(':')

    h,m,s=map(int,[h,m,s])
    # month,day,year=map(int,line.split()[0].split('/'))
    # h,m,s=map(int,line.split()[1].split(':'))
    # ampm=line.split()[2][:2]
    # # print(year,month,day,h,m,s)
    s0=calc_secs_explicit(h,m,s,day,year,ampm)
    return s0

def calc_vumat_pre(cnt,echo=False):
    lines=cnt.split('\n')
    i0=-1
    i1=0
    for i in range(len(lines)):
        if r'Begin Compiling Double Precision Abaqus/Explicit User Subroutines' in lines[i]:
            i0=i+1
        if r'End Abaqus/Explicit Packager' in lines[i]:
            i1=i-1

    # print('Estimating time for pre of EXPLICIT')
    # print('i0:',i0,' i1:',i1)

    # print(lines[i0])
    # print(lines[i1])

    try:
        s0=read_line_explicit_windows(lines[i0])
        s1=read_line_explicit_windows(lines[i1])
    except:
        s0=read_line_explicit_linux(lines[i0])
        s1=read_line_explicit_linux(lines[i1])

    return s1-s0

   # # dtime=time_delta(lines[i0],lines[i1])#/3600.

    # try:
    #     dtime=time_delta(lines[i0],lines[i1])#/3600.
    #     if echo:
    #         print('-------------------------------------------')
    #         for i in range(len(lines)):
    #             print(lines[i])
    #         print('-------------------------------------------')
    # except:
    #     print('Something must have gone wrong')
    #     print('Below is the actual lines from standard.exe')
    #     print('-------------------------------------------')
    #     if echo:
    #         for i in range(len(lines)):
    #             print(lines[i])
    #             print('-------------------------------------------')
    #     return np.nan
    # return dtime


def main(cnt,echo=False):
    """
    Arguments
    ---------
    cnt
    echo

    Returns
    -------
    Determine the total running time based on the
    message lines created from Abaqus.
    """
    import sys
    lines=cnt.split('\n')
    i0=-1
    i1=0

    beg_msg=['Run standard.exe','Run standard','Run explicit','Run package.exe']
    end_msg=['End Abaqus/Standard Analysis','End Abaqus/Explicit Analysis']

    for i, line in enumerate(lines):
        for j, phr in enumerate(beg_msg):
            if phr in line:
                i0=i-1
        for j, phr in enumerate(end_msg):
            if phr in line:
                i1=i-1

    ## in case failed to find the beginning
    if i0==-1:
        msg='** Error: could not determine i0, the'
        msg=f'{msg}\n   beginning of lines'
        sys.stderr.write(msg)
        sys.exit(-1)
        raise IOError(msg)

    ## in case failed to find the end
    if i1==0:
        finish=False
        k=0
        i=0
        while not(finish) and k<20:
            k=k+1
            i=i-1
            l=lines[i]
            if len(l)>6:
                dum=l[::-1][:6][::-1]
                if dum=='AM KST' or dum=='PM KST':
                    finish=True
                    i1=i
    ## in case still failed to find the end
    if i1==0:
        finish=False
        k=0
        i=0
        while not(finish) and k<20:
            k=k+1
            i=i-1
            l=lines[i]
            if len(l)>6:
                dum=l[::-1][:2][::-1]
                if dum=='AM' or dum=='PM':
                    finish=True
                    i1=i

    if i1==0:
        msg='** Error: could not determine i1, the'
        msg=f'{msg}\n   end of lines'
        sys.stderr.write(msg)
        sys.exit(-1)
        raise IOError(msg)

    # print(f'i0:  {i0} i1: {i1}')
    # print(f'lines[i0]: {lines[i0]}')
    dtime=time_delta(lines[i0],lines[i1])#/3600.
    try:
        dtime=time_delta(lines[i0],lines[i1])#/3600.
        if echo:
            print('-------------------------------------------')
            for i, line in enumerate(lines):
                print(line)
            print('-------------------------------------------')
    except:
        print('** Warning: something must have gone wrong!')
        print('Below is the actual lines from standard.exe')
        print('-------------------------------------------')
        if echo:
            for i in range(len(lines)):
                print(lines[i])
                print('-------------------------------------------')
        return np.nan
    return i0,i1,lines[i0:i1+1],dtime

def read_block(bl):
    lines=bl.split('\n')

    ## basic information
    line=lines[1]
    stp,inc,noel=map(int,line.split()[:3])
    dt=float(line.split()[9])

    ## time stamps
    t=lines[7] ## date and time

    rst=dict()
    rst.update(stp=stp,inc=inc,noel=noel,dt=dt,timestamp=t)
    return rst

def get_time(timestamp):
    import sys
    date,hms=timestamp.split()
    year,mon,day=map(int,(date[:4],date[4:6],date[6:8]))
    hrs,mins,secs=map(float,hms.split(':'))
    rst=dict(year=year,mon=mon,day=day,hrs=hrs,mins=mins,secs=secs)
    return rst

def get_dtime(timestamp1,timestamp2):
    """
    Calculate dtime between timestamp1 and timestamp2)
    """
    import sys
    rst1=get_time(timestamp1)
    rst2=get_time(timestamp2)

    if rst1['year']!=rst2['year']:
        msg='** Error: year difference was not considered'
        msg=f'{msg} in <calc_dtime.get_dtime>'
        sys.stderr.write(msg)
        sys.exit(-1)
    if rst1['mon']!=rst2['mon']:
        msg='** Error: month difference was not considered'
        msg=f'{msg} in <calc_dtime.get_dtime>'
        sys.stderr.write(msg)
        sys.exit(-1)

    dday=rst2['day']-rst1['day']
    dhrs=rst2['hrs']-rst1['hrs']
    dmin=rst2['mins']-rst1['mins']
    dsec=rst2['secs']-rst1['secs']

    dday=dday*24*3600
    dhrs=dhrs*3600
    dmin=dmin*60
    dtime=dday+dhrs+dmin+dsec
    return dtime

def find_wallclocktime_from_block(block,iopt=0):
    """
    Arguments
    ---------
    block
    iopt   (0: looks up backward)
           (1: looks up forward)
    """
    lines=block.split('\n')
    if iopt==1: lines=lines[::-1]
    timestamp=None
    for i, line in enumerate(lines):
        ## condition to define the time stamp
        if ':' in line and len(line.split())==2 and line[:2]=='20':
            timestamp=line
            break
    return timestamp

def find_fedt_from_block(block):
    lines=block.split('\n')
    items=lines[1].split()
    stp,inc,nel,npt,spt,de,edot,evm,svm,dt,dtdrv=items[:11]
    return float(dt),float(dtdrv)


def get_stats_timestamps(timestamps):
    """
    """
    dt_re_hist=np.zeros(len(timestamps)-1)
    for idt in range(len(timestamps)-1):
        t1,t2=timestamps[idt],timestamps[idt+1]
        dt_re_hist[idt]=get_dtime(t1,t2)
    # disp_stats_dt_hist(dt_re_hist)
    return dt_re_hist

def disp_stats_dt_hist(dt_hist):
    print(f'mean: {dt_hist.mean():.3f}')
    print(f'min : {dt_hist.min():.3f}')
    print(f'max : {dt_hist.max():.3f}')
    print(f'std : {dt_hist.std():.3f}')

def analysis(cnt):
    """
    Arguments
    ---------
    cnt
    """
    import sys

    i0,i1,lines,dtime=main(cnt)
    print(f'* Total FE run time: {dtime}')

    # dt in Real world.
    ## time (wall clock) elapsed for each attemp or each increment
    sep_msgs=['stpinc','** At the end   of current analysis increment']
    nbl=[]
    for imsg, sep_msg in enumerate(sep_msgs):
        blocks=cnt.split(sep_msg)
        timestamps=[]
        for ib, block in enumerate(blocks):
            if ib==len(blocks)-1: ## for the last one,
                timestamp=find_wallclocktime_from_block(block,iopt=1)
            else:
                timestamp=find_wallclocktime_from_block(block,iopt=0)
            if type(timestamp)!=type(None):
                timestamps.append(timestamp)
        print(f'\n** timestamp based on "{sep_msg[:6]}..."')
        print(f'len(timestamps):{len(timestamps)}')
        timestamps=np.array(timestamps,dtype='object')
        nbl.append(len(timestamps))
        dt_re_hist=get_stats_timestamps(timestamps)

    # print(f'nbl: {nbl}')
    n_total_attempts=nbl[0]
    n_total_incrs=nbl[1]
    print(f'Efficiency 1: {dtime/n_total_attempts:.3f}')
    print(f'Efficiency 2: {dtime/n_total_incrs:.3f}')
    print(f'Efficiency 3: {n_total_attempts/n_total_incrs:.3f}')

    # dt in FE world
    print(f'\n** timestamp based on "dt FE world"')
    blocks=cnt.split('stpinc')
    blocks=blocks[1:]
    fedt_hist=np.zeros(len(blocks))
    for ib, block in enumerate(blocks):
        fedt,fedtdrv=find_fedt_from_block(block)
        fedt_hist[ib]=fedt
    print(f'len(fedt_hist):{len(fedt_hist)}')
    disp_stats_dt_hist(fedt_hist)

def calc_time_trend(cnt,nol=30):
    length_checker = np.vectorize(len)
    ## summarize cnt
    lines_orig=cnt.split('\n')
    len_lines=length_checker(lines_orig)
    lines_orig=np.array(lines_orig,dtype=object)
    ind=len_lines==22 ## finding lines resulting from date and time
    lines=lines_orig[ind]
    D=[]

    print('---------------')
    print('a few last lines')
    for i in range(nol):
        print(lines_orig[-nol+i])
    print('---------------')

    if len(lines)<2:
        print('none')
        return
    else:
        for i in range(len(lines)-1):
            l0=lines[i]
            l1=lines[i+1]
            dtime=read_dt(l0,l1)
            D.append(dtime)

    D=np.array(D)
    print('---------------')
    print('mean: %6.3f [sec]'%D.mean())
    print('std:  %6.3f [sec]'%D.std())
    print('---------------')
    tt=read_dt(lines[0],lines[-1])
    print('dt:   %6.2f [min]   %6.1f [hr]   %7.1f [sec]'%(tt/60.,tt/3600.,tt))

    ## finding dtime(1)
    ind=len_lines==135 ## need to be estimated from the data
    lines=np.array(lines_orig,dtype=object)[ind]
    lines=lines_orig[ind]

    for i in range(len(lines)):
        l=lines[i]
        l=l.split('acctimeFE')[1]
        t1=float(l.split()[0])

        l=l.split('kinc')[1]
        kinc=float(l.split()[0])

    timespec=tt/t1
    timeinc=tt/kinc
    print('Calc time per 1 step time: %6.3f [hr] %5.0f [min] %6.0f [sec]'%(timespec/3600.,timespec/60.,timespec))
    print('Calc time per incr       : %6.3f [hr] %5.2f [min] %6.2f [sec]'%(timeinc/3600.,timeinc/60.,timeinc))

def read_dt(l0,l1):
    date1,secs1=read_time_from_line(l0)
    date2,secs2=read_time_from_line(l1)
    date_delta=date2-date1
    date_delta=date_delta.total_seconds()
    dtime = date_delta+secs2-secs1
    return dtime

def read_time_from_line(l='20210119  18:11:53.521'):
    from datetime import datetime,date

    time=l.split()[-1]
    h,m,s=time.split(':')
    h=int(h)
    m=int(m)
    secs=h*3600+m*60+float(s) #+ampm

    DD=l.split()[0]
    year,mon,day=int(DD[:4]),int(DD[4:6]),int(DD[6:8])
    date=date(year,mon,day)
    return date,secs

if __name__=='__main__':
    import argparse
    import sys
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', nargs='+', help='name of file',default='run.log')
    parser.add_argument(
        '-i', type=int,default=0,
        help='iopt (0: total time, 1: incremental time, 2: time / step analysis, 3: doing option 2 for multiple files)')
    args = parser.parse_args()

    if args.i==0 or args.i==1:
        for i, f in enumerate(args.f):
            with open(f,'r') as fo:
                cnt=fo.read()
            if args.i==0:
                try:
                    ## in case you are running vumat
                    dtime_prep=calc_vumat_pre(cnt,echo=False)
                    print('** dtime for preparation:',dtime_prep)
                except:
                    ## you are running umat
                    pass
                i0,i1,lines,dtime=main(cnt)
                print('** Total FE run time:',dtime)
            elif args.i==1:
                calc_time_trend(cnt)

    if args.i==2:
        for i, f in enumerate(args.f):
            with open(f,'r') as fo:
                cnt=fo.read()
            # i0,i1,lines,dtime=main(cnt)
            # print(f'dtime: {dtime}')
            rst=analysis(cnt)
