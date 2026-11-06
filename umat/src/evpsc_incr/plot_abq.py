import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from glob import glob
import os
def plot_vumat_vel(fn):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    dat=np.loadtxt(fn).T
    nd=len(dat[0])
    udot=np.zeros((3,3,nd))
    k=0
    for i in range(3):
        for j in range(3):
            #ax.plot(dat[0],dat[2+k])
            udot[i,j]=dat[2+k]
            k=k+1

    for i in range(3):
        for j in range(3):
            ax.plot(dat[0],udot[i,j],label=r'$u_{%i%i}$'%(i+1,j+1))

    ax.legend()
    fig.savefig('vel.pdf')

def plot_vumat_bm(path='/scratch/a1232a01/umat_scratches/vumat_benchmark'):
    fns=glob(os.path.join(path,'4_0008_100_*/umat/src/evpsc_incr/run.log'))
    fns.sort()

    fns_str=glob(os.path.join(path,'4_0008_100_*/umat/src/evpsc_incr/tmp_0001_0001.rst'))
    fns_str.sort()

    plt.ioff()
    fig=plt.figure(figsize=(12,6))
    ax1=fig.add_subplot(241)
    ax2=fig.add_subplot(242)
    ax3=fig.add_subplot(243)
    ax4=fig.add_subplot(244)
    ax5=fig.add_subplot(245)
    ax6=fig.add_subplot(246)

    x=[];y=[]
    for i in range(len(fns)):
        fn=fns[i]
        print('fn:',fn)
        try:
            msf=float(fn.split(path)[-1].split(os.sep)[1].split('_')[-1])
            cmd='python calc_dtime.py -f %s -i 0'%fn
            p=os.popen(cmd)
            cnt=p.read()
            dt=float(cnt.split('dtime:')[-1])
            print('dt:',dt)
            x.append(msf)
            y.append(dt)

            dat=np.loadtxt(fns_str[i],skiprows=1).T
            ax2.plot(dat[6],dat[12],label='%8.1e'%msf)
            ax3.plot(dat[6],dat[18],label='%8.1e'%msf)
            ax4.plot(dat[2][1:],np.diff(dat[2]),label='%8.1e'%msf)
        except:
            print('failed')
            pass

    ax1.plot(x,y,'kx')

    noels=[8,64,512,1000]

    symbols=['+','x','s','o']

    for j in range(len(noels)):
        fns_cpu=glob(os.path.join(path,'*_%4.4i_500_1e8/umat/src/evpsc_incr/run.log'%noels[j]))
        fns_cpu.sort()

        x=[];y=[]
        for i in range(len(fns_cpu)):
            fn=fns_cpu[i]
            print('fn_cpu:',fn)
            ncpu=float(fn.split(path)[-1].split(os.sep)[1].split('_')[0])

            try:
                cmd='python calc_dtime.py -f %s -i 0'%fn
                p=os.popen(cmd)
                cnt=p.read()
                dt=float(cnt.split('dtime:')[-1])
                print('dt:',dt)
                x.append(ncpu)
                y.append(dt)
            except:
                pass
        ax5.plot(x,y,symbols[j],label='NOEL: %i'%noels[j])
        ax6.plot(x,y,symbols[j],label='NOEL: %i'%noels[j])


    for ax in [ax5,ax6]:
        ax.set_xlabel('Number of CPU core units in use')
        ax.set_ylabel('Running time [s]')
        ax.legend()

    ax5.set_ylim(0.,)
    ax5.set_xlim(0.,)
    ax6.set_xscale('log',basex=2)
    ax6.set_yscale('log',basey=2)


    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(1e0,)

    ax1.set_xlabel('Mass Scale Factor')
    ax1.set_ylabel('Running time [sec]')

    ax2.legend()
    ax2.set_xlabel('Strain')
    ax2.set_ylabel('Stress [MPa]')

    ax3.set_xlabel('Step time [s]')
    ax3.set_ylabel(r'$\dot{\varepsilon}_{11}$ [/s]')
    ax4.set_ylabel('Time increment [s]')
    ax4.set_xlabel('Time [s]')

    fig.tight_layout()
    fig.savefig('vumat_bm.pdf')

class FErst:
    def __init__(self,fnrst,ntens=6):
        self.rawdata=np.loadtxt(fnrst,skiprows=1).T
        self.time=self.rawdata[2,:]
        self.xyz=self.rawdata[3:6,:]
        i0=6
        i1=i0+ntens
        self.LE=self.rawdata[i0:i1,:]
        i0=i1
        i1=i0+ntens
        self.S=self.rawdata[i0:i1,:]
        i0=i1
        i1=i0+ntens
        self.ER=self.rawdata[i0:i1,:]
        self.UVARM=self.rawdata[i1:,:]
        print(f'** self.UVARM.shape in <plot_abq.FErst>: {self.UVARM.shape}')
        if self.UVARM.shape[0]!=30:
            msg=f'** Error: the Number of UVARM is expected to be 30'
            sys.stderr.write(f'{msg}\n')
            sys.exit(-1)
            raise IOError(msg)

def plot_abq_rst(fn=None,ntens=6,fnsa=None,ihardlaw=0,fndd=None):
    """
    Arguments
    ---------
    fn   ('*.rst')
    ntens
    fnsa ('STR_STR.OUT')
    ihardlaw
    fndd
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from MP.mat.mech import FlowCurve as FC
    from plot import gen_fig
    import sys


    print(f'fn  : {fn}')
    print(f'fnsa: {fnsa}')
    # sys.exit(-1)

    plt.ioff()

    sa=FC()
    sa.get_model(fnsa)
    sa.time=np.cumsum(sa.tincrs)
    fe=FErst(fn,ntens)


    ncol=3
    if ihardlaw==0:
        nrow=1
    else:
        nrow=2

    fig,axes=gen_fig(ncols=ncol,nrows=nrow,
                     top=0.75,wspace=0.5,hspace=0.5,
                     bottom=0.3,right=0.8,colsize=3,rowsize=2.5)

    for i in range(3):
        kws=dict(label=r'$\dot{u}_{%i%i}$ SA'%(i+1,i+1),c=['r','g','b'][i],ls='-',lw=0.3,alpha=1.)
        fig.axes[0].plot(sa.time,sa.velgrads[i,i],**kws)
    ks=[0,4,8]
    for i, k in enumerate(ks):
        kws=dict(label=r'$\dot{u}_{%i%i}$ FE'%(i+1,i+1),c=['r','g','b'][i],ls='--',lw=3,alpha=0.3)
        fig.axes[0].plot(fe.time,fe.UVARM[k],**kws)

    ## von Mises stress & strain responses.
    fig.axes[1].plot(sa.epsilon_vm,sa.sigma_vm,'-',alpha=1,lw=1,label=r'$\Delta$EVPSC-SA')
    fig.axes[1].plot(fe.UVARM[9],fe.UVARM[10],'--',alpha=0.3,lw=3,mfc='None',label=r'$\Delta$EVPSC-FE')

    ## stress
    for i in range(3):
        kws=dict(c=['r','g','b','m'][i],ls='-',alpha=1.,lw=1)
        fig.axes[2].plot(sa.time,sa.sigma[i,i],label=r'$\sigma_{%i%i}, \Delta$EVPSC-SA'%(i+1,i+1),**kws)
        kws.update(ls='--',alpha=0.3,lw=3)
        fig.axes[2].plot(fe.time,fe.S[i],mfc='None',label=r'$\sigma_{%i%i}, \Delta$EVPSC-FE'%(i+1,i+1),**kws)

    ## if ihardlaw !=0
    if ihardlaw==20: ## dislocation density law of IB & CNT
        ax=fig.axes[3]
        nslmod=3 #(hardwired) Quite risky bet
        dat_dd=np.loadtxt(fndd,skiprows=1).T
        x=sa.time
        xfe=fe.time
        icol  =0   -1
        icoluv=11  -1
        ## forest dislocation per each slip mode

        for isl in range(nslmod):
            icol+=1
            icoluv+=1
            y=dat_dd[icol]
            kws=dict(c=['r','g','b','m'][isl],ls='-',lw=0.3,alpha=1)
            ax.plot(x,y,label=r'$\rho_{for}^{\alpha=%i}$, SA'%(isl+1),**kws)
            yfe=fe.UVARM[icoluv]
            kws.update(ls='--',lw=3,alpha=0.3)
            # print(yfe)
            ax.plot(fe.time,yfe,label=r'$\rho_{for}^{\alpha=%i}$, FE'%(isl+1),**kws)


        ax=fig.axes[4]
        axfe=fig.axes[5]
        # debris
        icol+=1
        icoluv+=1
        y=dat_dd[icol]
        yfe=fe.UVARM[icoluv]
        kws  =dict(label=r'$\rho_{deb}, \Delta$EVPSC-SA',c='r',ls='-',lw=0.3,alpha=1.0)
        kwsfe=dict(label=r'$\rho_{deb}, \Delta$EVPSC-FE',c='r',ls='--',lw=3.0,alpha=0.3)
        ax.plot(x,y,**kws)
        axfe.plot(xfe,yfe,**kwsfe)


        ## total
        icol+=1
        icoluv+=1
        y  =dat_dd[icol]
        yfe=fe.UVARM[icoluv]
        kws.update(label=r'$\rho_{tot}, \Delta$EVPSC-SA',c='g')
        kwsfe.update(label=r'$\rho_{tot}, \Delta$EVPSC-FE',c='r')
        ax.plot(x,y,**kws)
        axfe.plot(xfe,yfe,**kwsfe)

        ## debris per each slip mode
        for isl in range(nslmod):
            icol+=1
            icoluv+=1
            y=dat_dd[icol]
            yfe=fe.UVARM[icoluv]
            kws.update(  label=r'$\rho_{deb}^{\alpha=%i}, \Delta$EVPSC-SA'%(isl+1),c=['b','m','y'][isl])
            kwsfe.update(label=r'$\rho_{deb}^{\alpha=%i}, \Delta$EVPSC-FE'%(isl+1),c=['b','m','y'][isl])
            ax.plot(x,y,**kws)
            axfe.plot(xfe,yfe,**kwsfe)

        ## drho_recovery
        icol+=1
        icoluv+=1
        y=dat_dd[icol]
        yfe=fe.UVARM[icoluv]
        kws.update(  label=r'$\Delta \rho_{rec}, \Delta$EVPSC-SA',c='k')
        kwsfe.update(label=r'$\Delta \rho_{rec}, \Delta$EVPSC-FE',c='k')
        ax.plot(x,y,**kws)
        axfe.plot(xfe,yfe,**kwsfe)

    elif ihardlaw==23: ## dislocation density law of RGVB
        print('** Place holder for ihardlaw of 23 (RGVB)')
        pass



    ## -- deco begins
    fig.axes[0].set_xlabel('Time [s]')
    fig.axes[0].set_ylabel(r'$\dot{u}$ [/s]')
    fig.axes[1].set_xlabel(r'$\varepsilon^{VM}$')
    fig.axes[1].set_ylabel(r'$\sigma^{VM}$ [MPa]')
    fig.axes[2].set_xlabel('Time [s]')
    fig.axes[2].set_ylabel(r'$\sigma_{11},\sigma_{22},\sigma_{33}$ [MPa]')
    fig.axes[0].legend(loc='lower left',ncol=2,bbox_to_anchor=(0,1.05))
    fig.axes[2].legend(loc='lower left',ncol=1,bbox_to_anchor=(1.05,0))
    if ihardlaw==20:
        fig.axes[3].legend(loc='upper right',ncol=3,bbox_to_anchor=(1.2,-0.4),fontsize=7)
        fig.axes[3].set_yscale('log')
        fig.axes[4].legend(loc='upper left',ncol=3,bbox_to_anchor=(0,-0.4),fontsize=7)
        fig.axes[4].set_yscale('log')
        # fig.axes[5].legend(loc='upper left',ncol=3,bbox_to_anchor=(0,-0.4),fontsize=7)
        fig.axes[5].set_yscale('log')
    ## -- deco ends

    fnFigure_pdf=f"{fn.split('.')[0]}.pdf"
    fnFigure_png=f"{fn.split('.')[0]}.png"
    print(f'Figure saved to {fnFigure_pdf}')
    print(f'Figure saved to {fnFigure_png}')
    # fig.tight_layout()
    fig.savefig(fnFigure_pdf)
    fig.savefig(fnFigure_png)
    fig.savefig('tmp.pdf')
    fig.savefig('tmp.png')

if __name__=='__main__':
    import argparse
    ## Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',type=int,help='ioption, 0: plot using plot_abq_rst, 1: run <plot_vumat_vel> ',default=0)
    parser.add_argument('--ntens',type=int,help='ntens either 4 or 6')
    parser.add_argument('--ihardlaw',type=int,help='Hardening law',default=0)
    parser.add_argument('--fnRST',type=str,help='abq_rst file name',default='')
    parser.add_argument('--fnSA',type=str,help='dEVPSC stand-alone result',default='')
    parser.add_argument('--fnVEL',type=str,help='Velocity gradient',default='')
    parser.add_argument('--fndd',type=str,help='dislocation density file from Stand-alone calculation')

    args = parser.parse_args()

    if args.i==0:
        plot_abq_rst(fn=args.fnRST,ntens=args.ntens,fnsa=args.fnSA,ihardlaw=args.ihardlaw,fndd=args.fndd)
    elif args.i==1:
        plot_vumat_vel(args.fnVEL)
