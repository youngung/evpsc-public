# def rolling(edot,dt,eps,t0,t1,iopt=0):
#     """
#     rolling (plane-strain and no shear)

#     iopt=0 gives rolling for EVPSC model
#     iopt=1 gives rolling for VPSC model
#     """
#     import numpy as np
#     ## voigt component to cartesian i,j components
#     ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
#     iudot=np.ones((3,3),dtype='int')
#     iscau=np.zeros((3,3),dtype='int')
#     udot=np.zeros((3,3),dtype='float')
#     scauchy=np.zeros((3,3),dtype='float')
#     ctrlincr=abs(edot*dt)
#     nsteps=abs(int(eps/ctrlincr))
#     for v in range(6):
#         i,j=ijv[v]
#         if v+1==3:
#             iudot[i,j]=1
#             iscau[i,j]=0
#             udot[i,j]=-abs(edot)
#         elif v+1==2:
#             iudot[i,j]=1
#             iscau[i,j]=0
#             udot[i,j]=0.
#         else:
#             iscau[i,j]=1
#             iudot[i,j]=0
#     for v in np.arange(3,6):
#         i,j=ijv[v]
#         iudot[i,j]=1
#         iudot[j,i]=1
#         iscau[i,j]=0
#         iscau[j,i]=0
#     cnt=f'{nsteps} 3 {ctrlincr} {t0} {t1}'
#     cnt=f'{cnt}\n ** boundary conditions\n'
#     for i in range(3):
#         for j in range(3):
#             cnt=f'{cnt} %12.1i'%iudot[i,j]
#         cnt=f'{cnt}\n'
#     cnt=f'{cnt}\n'
#     for i in range(3):
#         for j in range(3):
#             cnt=f'{cnt} %12.3e'%udot[i,j]
#         cnt=f'{cnt}\n'
#     cnt=f'{cnt}\n'
#     for i in range(3):
#         for j in range(3):
#             if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
#             else   : cnt=f'{cnt} %12s'%''
#         cnt=f'{cnt}\n'
#     cnt=f'{cnt}\n'
#     for i in range(3):
#         for j in range(3):
#             if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
#             else   : cnt=f'{cnt} %12s'%''
#         cnt=f'{cnt}\n'
#     return cnt

def uniaxial_i(ictrl,edot,dt,eps,t0,t1):
    """
    uniaxial history file as function of
    <ictrl>, <edot>, <dt>, <eps>, <t0>, and <t1>
    """
    import numpy as np
    ## voigt component to cartesian i,j components
    ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
    iudot=np.ones((3,3),dtype='int')
    iscau=np.zeros((3,3),dtype='int')
    udot=np.zeros((3,3),dtype='float')
    scauchy=np.zeros((3,3),dtype='float')
    ctrlincr=abs(edot*dt)
    nsteps=int(eps/ctrlincr)
    for v in range(6): ## for all components
        i,j=ijv[v]
        if v+1==ictrl:
            iudot[i,j]=1
            iscau[i,j]=0
            udot[i,j]=edot
        else:
            iscau[i,j]=1
            iudot[i,j]=0

    for v in np.arange(3,6): ## for shear components
        i,j=ijv[v]
        if ictrl<=3:
            iudot[i,j]=1
            iudot[j,i]=0
            iscau[i,j]=1
            iscau[j,i]=1
        else:
            if v+1==ictrl:
                iudot[i,j]=1
                iudot[j,i]=1
                iscau[i,j]=0
                iscau[j,i]=0
            else:
                iudot[i,j]=1
                iudot[j,i]=0
                iscau[i,j]=1
                iscau[j,i]=1

    ## The below block is to make VPSC load_condition comparable,
    ## it should not affect the results of EVPSC calculations.
    if ictrl>=1 and ictrl <=3: ## uniaxial loading for normal component.
        for v in range(3):
            i,j=ijv[v]
            if v+1!=ictrl:
                udot[i,j]=-edot/2.


    cnt=f'{nsteps} {ictrl} {ctrlincr} {t0} {t1}'
    cnt=f'{cnt}\n ** boundary conditions\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.1i'%iudot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.3e'%udot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    return cnt



def bulge(edot,dt,eps,t0,t1):
    """
    uniaxial history file as function of
    <ictrl>, <edot>, <dt>, <eps>, <t0>, and <t1>
    """
    import numpy as np
    ictrl=1

    ## voigt component to cartesian i,j components
    ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
    iudot=np.ones((3,3),dtype='int')
    iscau=np.zeros((3,3),dtype='int')
    udot=np.zeros((3,3),dtype='float')
    scauchy=np.zeros((3,3),dtype='float')
    ctrlincr=abs(edot*dt)
    nsteps=int(eps/ctrlincr)
    for v in range(6):
        i,j=ijv[v]
        if v+1==1 or v+1==2:
            iudot[i,j]=1
            iscau[i,j]=0
            udot[i,j]=edot
        else:
            iscau[i,j]=1
            iudot[i,j]=0
    for v in np.arange(3,6):
        i,j=ijv[v]
        iudot[i,j]=1
        iudot[j,i]=0
        iscau[i,j]=1
        iscau[j,i]=1
    cnt=f'{nsteps} {ictrl} {ctrlincr} {t0} {t1}'
    cnt=f'{cnt}\n ** boundary conditions\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.1i'%iudot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.3e'%udot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    return cnt


def rolling(edot,dt,eps,t0,t1,iopt=0):
    """
    rolling (plane-strain and no shear)

    iopt=0 gives rolling for EVPSC model
    iopt=1 gives rolling for VPSC model
    """
    import numpy as np
    ictrl=3

    ## voigt component to cartesian i,j components
    ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
    iudot=np.ones((3,3),dtype='int')
    iscau=np.zeros((3,3),dtype='int')
    udot=np.zeros((3,3),dtype='float')
    scauchy=np.zeros((3,3),dtype='float')
    ctrlincr=abs(edot*dt)
    nsteps=int(eps/ctrlincr)
    for v in range(6):
        i,j=ijv[v]
        if v+1==3: ## thick
            iudot[i,j]=1
            iscau[i,j]=0
            udot[i,j]=edot
        elif v+1==2 and iopt==0: ## transverse
            iudot[i,j]=0
            iscau[i,j]=1
            udot[i,j]=0.
        elif v+1==2 and iopt==1: ## transverse
            iudot[i,j]=1
            iscau[i,j]=0
            udot[i,j]=0.
        elif v+1==1: ## rolling
            iudot[i,j]=1
            iscau[i,j]=0
            udot[i,j]=-edot
        else:
            iudot[i,j]=0
            iscau[i,j]=1
            udot[i,j]=0.

    cnt=f'{nsteps} {ictrl} {ctrlincr} {t0} {t1}'
    cnt=f'{cnt}\n ** boundary conditions\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.1i'%iudot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.3e'%udot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    return cnt

def psrd(edot,dt,eps,t0,t1):
    """
    plane strain RD
    <edot>, <dt>, <eps>, <t0>, and <t1>
    """
    import numpy as np
    ictrl=1

    ## voigt component to cartesian i,j components
    ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
    iudot=np.ones((3,3),dtype='int')
    iscau=np.zeros((3,3),dtype='int')
    udot=np.zeros((3,3),dtype='float')
    scauchy=np.zeros((3,3),dtype='float')
    ctrlincr=abs(edot*dt)
    nsteps=int(eps/ctrlincr)

    iudot[:,:]=0
    iscau[:,:]=0
    iudot[0,0]=1
    iudot[1,1]=1
    iudot[2,2]=0
    iudot[0,1]=0
    iudot[1,0]=1
    iudot[0,2]=0
    iudot[2,0]=1
    iudot[1,2]=0
    iudot[2,1]=1
    udot[0,0]=edot

    iscau[:,:]=0
    for v in range(6):
        i,j=ijv[v]
        if iudot[i,j]*iudot[i,j]==0:
            iscau[i,j]=1

    cnt=f'{nsteps} {ictrl} {ctrlincr} {t0} {t1}'
    cnt=f'{cnt}\n ** boundary conditions\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.1i'%iudot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.3e'%udot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    return cnt


def temp(time,dt,t0,t1):
    import numpy as np
    ## voigt component to cartesian i,j components
    ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
    ictrl=7
    iudot=np.ones((3,3),dtype='int')
    iscau=np.zeros((3,3),dtype='int')
    udot=np.zeros((3,3),dtype='float')
    scauchy=np.zeros((3,3),dtype='float')
#    scauchy[0,0]=1e-5 #np.sign(t1-t0)*(1e-5) ##minimal stress increment (to avoid num. error)
    scauchy[0,0]=np.sign(t1-t0)*(1e0) ##minimal stress increment (to avoid num. error)
    iudot[:,:]=0
    iudot[0,1]=1
    iudot[0,2]=1
    iudot[1,2]=1
    iscau[:,:]=1
    nsteps=int(time/dt)
    cnt=f'{nsteps} {ictrl} {dt} {t0} {t1}'
    cnt=f'{cnt}\n ** boundary conditions\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.1i'%iudot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.3e'%udot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'

    return cnt


def hold_eps(t0,t1,dt,time):
    """
    Hold the strain and let the stress of material relaxes.

    Arguments
    ---------
    t0
    t1
    dt
    time

    Returns
    -------
    cnt
    """
    import numpy as np
    ## voigt component to cartesian i,j components
    ijv=[[0,0],[1,1],[2,2],[1,2],[0,2],[0,1]]
    iudot=np.ones((3,3),dtype='int')
    iscau=np.zeros((3,3),dtype='int')
    udot=np.zeros((3,3),dtype='float')
    scauchy=np.zeros((3,3),dtype='float')
    nsteps=int(time/dt)
    ictrl=7

    iscau[:,:]=0
    iudot[:,:]=1

    cnt=f'{nsteps} {ictrl} {dt} {t0} {t1}'
    cnt=f'{cnt}\n ** boundary conditions\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.1i'%iudot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            cnt=f'{cnt} %12.3e'%udot[i,j]
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.1i'%iscau[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    cnt=f'{cnt}\n'
    for i in range(3):
        for j in range(3):
            if j>=i: cnt=f'{cnt} %12.3e'%scauchy[i,j]
            else   : cnt=f'{cnt} %12s'%''
        cnt=f'{cnt}\n'
    return cnt




if __name__=='__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ictrl',type=int,default=1)
    parser.add_argument(
        '--edot',type=float,default=1e-3)
    parser.add_argument(
        '--dt',type=float,default=1.0)
    parser.add_argument(
        '--eps',type=float,default=0.01,help='total strain / total time')
    parser.add_argument(
        '--t0',type=float,default=298.,help='Initial temperature')
    parser.add_argument(
        '--t1',type=float,default=298.,help='final temperature')
    parser.add_argument(
        '--iopt',type=int,default=0,
        help='0: uniaxial mechanical;  1: thermal (no strain); 2: hold (strain); 3: rolling; 4: bulge; 5: PSRD; 6: hold (stress); -3')
    parser.add_argument(
        '--fnout',type=str,default='',help='output file name')
    args=parser.parse_args()
    if args.iopt==0:
        hist_cnt=uniaxial_i(
            ictrl=args.ictrl,edot=args.edot,
            dt=args.dt,eps=args.eps,t0=args.t0,t1=args.t1)
    elif args.iopt==1:
        hist_cnt=temp(time=args.eps,dt=args.dt,t0=args.t0,t1=args.t1)
    elif args.iopt==2:
        hist_cnt=hold_eps(t0=args.t0,t1=args.t1,dt=args.dt,time=args.eps)
    elif args.iopt==3:
        hist_cnt=rolling(edot=args.edot,
                         dt=args.dt,eps=args.eps,t0=args.t0,t1=args.t1,iopt=0)
    elif args.iopt==-3:
        hist_cnt=rolling(edot=args.edot,
                         dt=args.dt,eps=args.eps,t0=args.t0,t1=args.t1,iopt=1)
    elif args.iopt==4:
        hist_cnt=bulge(
                       edot=args.edot,
                       dt=args.dt,eps=args.eps,t0=args.t0,t1=args.t1)
    elif args.iopt==5:
        hist_cnt=psrd(
            edot=args.edot,dt=args.dt,eps=args.eps,t0=args.t0,t1=args.t1)
    else:
        raise IOError('** Error: Unexpected iopt given')

    # print(hist_cnt, file=sys.stdout)
    if args.fnout!='':
        with open(args.fnout,'w') as fo:
            fo.write(hist_cnt)
