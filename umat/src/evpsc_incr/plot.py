## Plotting the simulated results from delta EVPSC-UMAT

from TX import upf
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from MP import read_blocks as rb
from MP.mat.mech import FlowCurve as FC

def gen_fig(nrows=3,ncols=3,colsize=2.5,rowsize=2.,**kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    gs=GridSpec(nrows=nrows,ncols=ncols,**kwargs)
    fig=plt.figure(figsize=(colsize*ncols,rowsize*nrows))
    axes=np.empty((nrows,ncols),dtype='object')
    for i in range(nrows):
        for j in range(ncols):
            ax=fig.add_subplot(gs[i,j])
            axes[i,j]=ax
    return fig,axes

def compare(fnLOG='dev0001_0001.log',fnSA='STR_STR.OUT',fnFE=None,ntens=6):
    """
    Arguments
    ---------
    fnLOG
    fnSA
    fnFE
    ntens
    """
    dat_LOG=rb.rb(fnLOG,1)
    print('fnLOG:',fnLOG)
    print('dat_LOG.shape:',dat_LOG.shape)
    fc_SA=FC()
    fc_SA.get_model(fnSA)
    fc_SA.time=np.cumsum(fc_SA.tincrs)

    fig=plt.figure(figsize=(16,16))
    ax1=fig.add_subplot(5,5,1)
    ax2=fig.add_subplot(5,5,2)
    ax3=fig.add_subplot(5,5,3)
    ax4=fig.add_subplot(5,5,4)
    ax5=fig.add_subplot(5,5,5)
    ax6=fig.add_subplot(5,5,6)
    ax7=fig.add_subplot(5,5,7)
    ax8=fig.add_subplot(5,5,8)
    ax9=fig.add_subplot(5,5,9)
    ax10=fig.add_subplot(5,5,10)
    ax11=fig.add_subplot(5,5,11)
    ax12=fig.add_subplot(5,5,12)
    ax13=fig.add_subplot(5,5,13)
    ax19=fig.add_subplot(5,5,19)
    ax14=fig.add_subplot(5,5,14)
    ax15=fig.add_subplot(5,5,15)



    if fnFE!='None':
        ax16=fig.add_subplot(5,5,16)
        ax17=fig.add_subplot(5,5,17)
        ax20=fig.add_subplot(5,5,20)
        ax21=fig.add_subplot(5,5,21)
        ax22=fig.add_subplot(5,5,22)



    ## von Mises stress-strain curve (EVM vs. SVM)
    ax1.plot(dat_LOG[30],dat_LOG[31],'-',label='LOG')
    ax1.plot(fc_SA.epsilon_vm,fc_SA.sigma_vm,'--',mfc='None',label='SA')

    ## S11-E11 stran11 vs. stress11
    ax2.plot(dat_LOG[9],dat_LOG[15],'-',label='LOG')
    ax2.plot(fc_SA.epsilon[0,0],fc_SA.sigma[0,0],'--',mfc='None',label='SA')
    ## S22-E22 stran22 vs. stress22
    ax10.plot(dat_LOG[10],dat_LOG[16],'-',label='LOG')
    ax10.plot(fc_SA.epsilon[1,1],fc_SA.sigma[1,1],'--',mfc='None',label='SA')

    ## S12-E12 stran12 vs. stress12
    ax3.plot(dat_LOG[11],dat_LOG[17],'-',label='LOG')
    ax3.plot(fc_SA.epsilon[0,1],fc_SA.sigma[0,1],'--',mfc='None',label='SA')

    # strain rate (LOG) time(2) vs. sr11, sr22, sr12
    ax4.plot(dat_LOG[7],dat_LOG[32],'r-',label=r'$\bar{\dot{\varepsilon}}_{11}$ (LOG)')
    ax4.plot(dat_LOG[7],dat_LOG[33],'g-',label=r'$\bar{\dot{\varepsilon}}_{22}$ (LOG)')
    ax4.plot(dat_LOG[7],dat_LOG[34],'b-',label=r'$\bar{\dot{\varepsilon}}_{12}$ (LOG)')

    # strain rate (SA)
    sr12=(fc_SA.velgrads[0,1]+fc_SA.velgrads[1,0])/2.
    sr22=fc_SA.velgrads[1,1]
    sr11=fc_SA.velgrads[0,0]
    ax4.plot(fc_SA.time[::10],sr11[::10],'ro',mfc='None',label=r'$\bar{\dot{\varepsilon}}_{11}$ (SA)')
    ax4.plot(fc_SA.time[::10],sr22[::10],'go',mfc='None',label=r'$\bar{\dot{\varepsilon}}_{22}$ (SA)')
    ax4.plot(fc_SA.time[::10],sr12[::10],'bo',mfc='None',label=r'$\bar{\dot{\varepsilon}}_{12}$ (SA)')

    # stress components (LOG), time(2) vs. stress11, stress22, stress 12
    ax5.plot(dat_LOG[7],dat_LOG[15],'r-',label=r'$\bar{\sigma}_{11}$ (LOG)')
    ax5.plot(dat_LOG[7],dat_LOG[16],'g-',label=r'$\bar{\sigma}_{22}$ (LOG)')
    ax5.plot(dat_LOG[7],dat_LOG[17],'b-',label=r'$\bar{\sigma}_{12}$ (LOG)')

    # stress components (SA)
    ax5.plot(fc_SA.time[::10],fc_SA.sigma[0,0][::10],'ro',mfc='None',label=r'$\bar{\sigma}_{11}$ (SA)')
    ax5.plot(fc_SA.time[::10],fc_SA.sigma[1,1][::10],'go',mfc='None',label=r'$\bar{\sigma}_{22}$ (SA)')
    ax5.plot(fc_SA.time[::10],fc_SA.sigma[0,1][::10],'bo',mfc='None',label=r'$\bar{\sigma}_{12}$ (SA)')

    # Modulus time(2) vs. J11, J22, J12
    ax15.plot(dat_LOG[7],dat_LOG[18]*1e-3,'r--',label=r'$\frac{d\Delta\bar{\sigma}_{1}}{d\Delta\bar{\epsilon}_{1}}$')
    ax15.plot(dat_LOG[7],dat_LOG[19]*1e-3,'g--',label=r'$\frac{d\Delta\bar{\sigma}_{2}}{d\Delta\bar{\epsilon}_{2}}$')
    ax15.plot(dat_LOG[7],dat_LOG[20]*1e-3,'b--',label=r'$\frac{d\Delta\bar{\sigma}_{1}}{d\Delta\bar{\epsilon}_{2}}$')

    ## SA analysis
    ax6.plot(fc_SA.time,fc_SA.sigma[0,0],label=r'$\bar{\sigma}_{11}$ (SA)')
    ax6.plot(fc_SA.time,fc_SA.sigma[1,1],label=r'$\bar{\sigma}_{22}$ (SA)')
    ax6.plot(fc_SA.time,fc_SA.sigma[2,2],label=r'$\bar{\sigma}_{33}$ (SA)')
    ax6.plot(fc_SA.time,fc_SA.sigma[0,1],label=r'$\bar{\sigma}_{12}$ (SA)')
    ax6.plot(fc_SA.time,fc_SA.sigma[0,2],label=r'$\bar{\sigma}_{13}$ (SA)')
    ax6.plot(fc_SA.time,fc_SA.sigma[1,2],label=r'$\bar{\sigma}_{23}$ (SA)')



    if ntens==6:
        ivo=[[0,0],[1,1],[2,2],[0,1],[0,2],[1,2]]
    elif ntens==3:
        ivo=[[0,0],[1,1],[0,1]]

    ## LOG and SA comparison on
    for i in range(ntens):
        k,l=ivo[i][0],ivo[i][1]
        # print('k,l:',k,l)
        # time(2) vs. eps_e0_v1, eps_e0_v2, eps_e0_v3, eps_e0_v4, eps_e0_v5, eps_e0_v6
        if ntens==6: l_,=ax12.plot(dat_LOG[7],dat_LOG[36+i],'-',label=r'$\bar{\epsilon}_{%i}$ (LOG)'%(i+1))
        if ntens==3: l_,=ax12.plot(dat_LOG[7],dat_LOG[35+i],'-',label=r'$\bar{\epsilon}_{%i}$ (LOG)'%(i+1))
        ax7.plot(fc_SA.time,fc_SA.epsilon[k,l],'-',mfc='None',label=r'$\bar{\epsilon}_{%i%i}$ (SA)'%(k+1,l+1),mec=l_.get_color())
        # time(2) vs. sig_e0_v1, sig_e0_v2, sig_e0_v3, sig_e0_v4, sig_e0_v5, sig_e0_v6
        if ntens==6:ax11.plot(dat_LOG[7],dat_LOG[42+i],'-',label=r'$\bar{\sigma}^{e^{(0)}}_{%i}$ (LOG)'%(i+1))
        if ntens==3:ax11.plot(dat_LOG[7],dat_LOG[41+i],'-',label=r'$\bar{\sigma}^{e^{(0)}}_{%i}$ (LOG)'%(i+1))

    if ntens==6:
        k=0
        for i in range(3):
            for j in range(3):
                ax13.plot(dat_LOG[7],dat_LOG[87+k],'-',label=r'$\bar{L}^{%i%i}$ (LOG)'%(i+1,j+1))
                k=k+1
        pass
    # in the case of ntens==6
    #dat_LOG[42]: stress
    #dat_LOG[60]: DFGRD0
    #dat_LOG[69]: DFGRD1
    #dat_LOG[78]: rot_acum ! e0 to e(n+1)
    #dat_LOG[87]: udot


    ## velgrads
    fc_SA.dsim=np.zeros(fc_SA.velgrads.shape)
    fc_SA.omega=np.zeros(fc_SA.velgrads.shape)
    fc_SA.dsim[:,:,:]=0.
    fc_SA.omega[:,:,:]=0.
    for i in range(3):
        for j in range(3):
            ax8.plot(fc_SA.time,fc_SA.velgrads[i,j],label=r'$\bar{L}_{%1i%1i}$ (SA)'%(i+1,j+1))
            fc_SA.dsim[i,j]=0.5*(fc_SA.velgrads[i,j]+fc_SA.velgrads[j,i])
            fc_SA.omega[i,j]=0.5*(fc_SA.velgrads[i,j]-fc_SA.velgrads[j,i])


    ## strain rate (SA)
    for i in range(ntens):
        k,l=ivo[i][0],ivo[i][1]
        ax19.plot(fc_SA.time,fc_SA.dsim[k,l],label=r'$\bar{\dot{\varepsilon}}_{%1i%1i}$ (SA)'%(k+1,l+1))

    ## spin (SA)
    for k in range(3):
        for l in range(3):
            ax14.plot(fc_SA.time,fc_SA.omega[k,l],label=r'$\bar{\dot{\omega}}_{%1i%1i}$ (SA)'%(k+1,l+1))


    ## r-value?
    ax9.plot(fc_SA.time,fc_SA.velgrads[1,1]/fc_SA.velgrads[2,2],
             label=r'$\bar{L}_{22}/\bar{L}_{33}$ (SA)')
    ax9.plot(fc_SA.time,fc_SA.velgrads[0,0]/fc_SA.velgrads[2,2],
             label=r'$\bar{L}_{11}/\bar{L}_{33}$ (SA)')


    if fnFE!='None':
        dat=np.loadtxt(fnFE,skiprows=1).T

        if ntens==6: nt = 6
        if ntens==3: nt = 4

        time_FE=dat[2]
        LE=dat[6:6+nt]
        S =dat[6+nt:6+nt*2]
        ER=dat[6+nt*2:6+nt*3]
        for i in range(nt):
            ax17.plot(time_FE,LE[i],label=r'$LE%i$'%(i+1))
        for i in range(nt):
            ax16.plot(time_FE,S[i],label=r'$S%i$'%(i+1))
        for i in range(nt):
            ax20.plot(time_FE[1:],ER[i][1:],label=r'$ER%i$'%(i+1))



        if ntens==6:
            #dat_LOG[78]: rot_acum
            S_33=np.zeros((3,3,dat.shape[-1]))
            LE_33=np.zeros((3,3,dat.shape[-1]))

            for i in range(6):
                k,l=ivo[i][0],ivo[i][1]
                S_33[k,l]=S[i]
                S_33[l,k]=S[i]
                LE_33[k,l]=LE[i]
                LE_33[l,k]=LE[i]


            R_33=np.zeros((3,3,dat_LOG.shape[-1]))
            k=0
            for i in range(3):
                for j in range(3):
                    R_33[i,j]=dat_LOG[78+k]
                    k=k+1

            print('R_33.shape:',R_33.shape)
            steps=[]
            incrs=[]
            R=[]
            for i in range(dat_LOG.shape[-1]):
                stp_=dat_LOG[2,i]
                inc_=dat_LOG[3,i]
                if inc_ in incrs and stp_ in steps:
                    incrs[-1]=inc_
                    steps[-1]=stp_
                    R[-1]=R_33[:,:,-1]
                else:
                    steps.append(stp_)
                    incrs.append(inc_)
                    R.append(R_33[:,:,i])

            R_old=np.array(R)
            R_new=np.zeros((3,3,R_old.shape[0]))
            for i in range(len(R_old)):
                R_new[:,:,i]=R_old[i,:,:]
            R=R_new[::]
            print('R.shape:',R.shape)
            steps=np.array(steps)
            print('steps.shape:',steps.shape)
            print('dat.shape:',dat.shape)


            S_=np.zeros(S_33.shape)
            LE_=np.zeros(LE_33.shape)
            for istp in range(dat.shape[-1]):
                try:
                    r=R[:,:,istp+1]
                    s=S_33[:,:,istp]
                    le=LE_33[:,:,istp]
                    S_[:,:,istp]=np.dot(r.T,np.dot(s,r))
                    LE_[:,:,istp]=np.dot(r.T,np.dot(le,r))
                except: pass

            # for i in range(6):
            #     k,l=ivo[i][0],ivo[i][1]
            #     ax21.plot(S_[k,l])
            #     ax22.plot(LE_[k,l])


        elif ntens==3:
            S_33=np.zeros((3,3,dat.shape[-1]))
            LE_33=np.zeros((3,3,dat.shape[-1]))

            for i in range(3):
                pass









    ## decorations
    ax1.set_xlabel(r'$\bar{\varepsilon}^{VM}$')
    ax2.set_xlabel(r'$\bar{\varepsilon}_{11}$ ')
    ax10.set_xlabel(r'$\bar{\varepsilon}_{22}$ ')
    ax3.set_xlabel(r'$\bar{\varepsilon}_{12}$')
    ax4.set_xlabel(r'time [sec]')
    ax5.set_xlabel(r'time [sec]')
    ax15.set_xlabel(r'time [sec]')
    ax1.set_ylabel(r'$\bar{\sigma}^{VM}$ [MPa]')
    ax2.set_ylabel(r'$\bar{\sigma}_{11}$ [MPa]')
    ax10.set_ylabel(r'$\bar{\sigma}_{22}$ [MPa]')
    ax3.set_ylabel(r'$\bar{\sigma}_{12}$ [MPa]')
    ax4.set_ylabel(r'Strain rate $\bar{\dot{\varepsilon}}_{ij} \ [s^{-1}]$')
    ax5.set_ylabel(r'Stress $\bar{\sigma}_{ij}$ [MPa]')
    ax15.set_ylabel(r'Jacobian [GPa]')
    ax6.set_ylabel(r'Stress $\bar{\sigma}_{ij}$ [MPa]')
    ax6.set_xlabel(r'time [sec]')
    ax11.set_ylabel(r'Stress $\bar{\sigma}_{ij}$ [MPa]')
    ax11.set_xlabel(r'time [sec]')
    ax7.set_ylabel(r'Strain $\bar{\epsilon}_{ij}$')
    ax7.set_xlabel(r'time [sec]')
    ax12.set_ylabel(r'Strain $\bar{\epsilon}_{ij}$')
    ax12.set_xlabel(r'time [sec]')
    ax8.set_ylabel(r'Velocity Gradient $\bar{L}_{ij} \ [s^{-1}]$')
    ax8.set_xlabel(r'time [sec]')
    ax13.set_ylabel(r'Velocity Gradient $\bar{L}_{ij} \ [s^{-1}]$')
    ax13.set_xlabel(r'time [sec]')
    ax9.set_xlabel(r'time [sec]')
    ax19.set_ylabel(r'Strain rate $\bar{\dot{\varepsilon}}_{ij} \ [s^{-1}]$')
    ax19.set_xlabel(r'time [sec]')
    ax14.set_ylabel(r'Spin $\bar{\dot{\omega}}_{ij} \ [s^{-1}]$')
    ax14.set_xlabel(r'time [sec]')

    if fnFE!='None':

        ax16.set_ylabel(r'S')
        ax16.set_xlabel(r'time [sec]')
        ax17.set_ylabel(r'LE')
        ax17.set_xlabel(r'time [sec]')
        ax20.set_ylabel(r'ER')
        ax20.set_xlabel(r'time [sec]')

    ax15.set_ylim(0,)
    ax9.set_ylim(0,)

    for ax in fig.axes[:14]:
        ax.legend(fontsize=6)

    if fnFE!='None':
        ax16.legend(fontsize=6)
        ax17.legend(fontsize=6)
        ax20.legend(fontsize=6)

    fnrst='figs/fig_compare.png'
    fig.tight_layout()
    fig.savefig(fnrst,dpi=200,bbonx_to_inches='tight')
    print(f"** {fnrst} has been created.")

    if True:
        # texture
        fns=glob('TEX_PH?.OUT')
        fns.sort()
        fns=np.array(fns,dtype='object')
        fnrst=np.empty(fns.shape,dtype='object')

        for i in range(len(fns)):
            grains,ngr=rb.read_tx(fns[i])
            px=grains[-1].T
            # px=grains[j].T
            mypf=upf.polefigure(grains=px,csym='cubic',cdim=[1,1,1],cang=[90,90,90])
            fig=mypf.pf_new(poles=[[1,0,0],[1,1,0],[1,1,1]])
            fnrst=f'pf_ph{i}.png'
            fig.savefig(fnrst)
            print(f'** {fnrst} has been created.')

def plot_texture(fn='TEX_PH1.OUT',poles=[[0,0,0,2]],**kwargs):
    """
    Arguments
    ---------
    file name of texture
    **kwargs passed to upf.polefigure class
    """
    from TX import upf
    pxs,ngrs=rb.read_tx(fn)

    #mypf=upf.polefigure(filename='TEX_PH1.OUT',**kwargs)

    mypf=upf.polefigure(pxs[0].T,**kwargs)
    fig=mypf.pf_new(poles=poles)
    fig.savefig('figs/polefigure_0.png')


    mypf=upf.polefigure(pxs[-1].T,**kwargs)
    fig=mypf.pf_new(poles=poles)
    fig.savefig('figs/polefigure_1.png')


def main(dst):#fn=None):
    import os
    from pathlib import Path

    if type(dst)==type(None):
        raise IOError('Unexected argument <fn> given')

    fn=os.path.join(dst,'dev0001_0001.log')
    dat=rb.main(fn)

    fig,axes=gen_fig(nrows=3,ncols=3,colsize=3,rowsize=3)

    ax1,ax2,ax3=axes[0,:3]
    ax4,ax5,ax6=axes[1,:3]

    ax1.plot(dat[30],dat[31],'-x')
    ax2.plot(dat[9], dat[15],'-x')
    ax3.plot(dat[11],dat[17],'-x')
    # strain rate
    ax4.plot(dat[7],dat[32],'-x',label=r'$\bar{\dot{\varepsilon}}_{11}$')
    ax4.plot(dat[7],dat[33],'-x',label=r'$\bar{\dot{\varepsilon}}_{22}$')
    ax4.plot(dat[7],dat[34],'-x',label=r'$\bar{\dot{\varepsilon}}_{12}$')
    # stress components
    ax5.plot(dat[7],dat[15],'-x',label=r'$\bar{\sigma}_{11}$')
    ax5.plot(dat[7],dat[16],'-x',label=r'$\bar{\sigma}_{22}$')
    ax5.plot(dat[7],dat[17],'-x',label=r'$\bar{\sigma}_{12}$')

    ax1.set_xlabel(r'$\bar{\varepsilon}^{VM}$')
    ax1.set_ylim(0,); ax2.set_ylim(0,)
    ax2.set_xlabel(r'$\bar{\varepsilon}_{11}$ ')
    ax3.set_xlabel(r'$\bar{\varepsilon}_{12}$')
    ax4.set_xlabel(r'time [sec]')
    ax5.set_xlabel(r'time [sec]')
    ax1.set_ylabel(r'$\bar{\sigma}^{VM}$ [MPa]')
    ax2.set_ylabel(r'$\bar{\sigma}_{11}$ [MPa]')
    ax3.set_ylabel(r'$\bar{\sigma}_{12}$ [MPa]')
    ax4.set_ylabel(r'Strain rate $\bar{\dot{\varepsilon}}_{ij} \ [s^{-1}]$')
    ax5.set_ylabel(r'Stress $\bar{\sigma}_{ij}$ [MPa]')

    ax4.legend(); ax5.legend()



    Path('figs').mkdir(exist_ok=True)
    fnrst='figs/fig.png'
    fig.tight_layout()
    fig.savefig(fnrst,dpi=200)#,bbox_to_inches='tight')
    print(f"** {fnrst} has been created.")

    if False:
    #if True:
        fns=glob(os.path.join(dst,'TEX_PH1_EL0000001  _KT0000000  .OUT'))
        print(f'fns: {fns}')
        fns.sort()
        fns=np.array(fns,dtype='object')
        fnrst=np.empty(fns.shape,dtype='object')
        for i in range(len(fns[:1])):
            print('fns[i]:',fns[i])
            grains,ngr=rb.read_tx(fns[i])
            px_last=grains[-1].T
            mypf_last=upf.polefigure(grains=px_last,csym='cubic',cdim=[1,1,1],cang=[90,90,90])
            px_first=grains[0].T
            mypf_first=upf.polefigure(grains=px_first,csym='cubic',cdim=[1,1,1],cang=[90,90,90])

            fig_last=mypf_last.pf_new(poles=[[1,0,0],[1,1,0],[1,1,1]])
            fnrst_last=f"figs/{fns[i].split('.OUT')[0]}_last.png"
            fig_last.savefig(fnrst_last)

            fig_first=mypf_first.pf_new(poles=[[1,0,0],[1,1,0],[1,1,1]])
            fnrst_first=f"figs/{fns[i].split('.OUT')[0]}_first.png"
            fig_first.savefig(fnrst_first)

            print(f'** {fnrst_first} and {fnrst_last} have been created.')

    if False:
        # texture
        fns=glob('TEX_PH?.OUT')
        fns.sort()
        fns=np.array(fns,dtype='object')
        fnrst=np.empty(fns.shape,dtype='object')

        for i in range(len(fns)):
            grains,ngr=rb.read_tx(fns[i])
            for j in range(len(grains)):
                px=grains[j].T

                mypf=upf.polefigure(grains=px,csym='cubic',cdim=[1,1,1],cang=[90,90,90])
                fig=mypf.pf_new(poles=[[1,0,0],[1,1,0],[1,1,1]])
                fnrst=f'figs/pf_{i}_{j}.png'
                fig.savefig(fnrst)

                print(f'** {fnrst} has been created.')

if __name__=='__main__':
    import argparse
    ## Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--i',type=int,help='option (0: main;  1: compare)')
    parser.add_argument('--dst',type=str,help='<dst>')
#    parser.add_argument('--fnLOG',type=str,help='vel. log file from FE EVPSC',default='dev0001_0001.log')
    parser.add_argument('--fnSA',type=str,help='STR_STR.OUT file from SA EVPSC',default='STR_STR.OUT')
    parser.add_argument('--fnFE',type=str,help='File generated from <abq_rst.py>',default='None')
    parser.add_argument('--ntens',type=int,help='NTENS',default=6)

    args = parser.parse_args()

    if args.i==0:
        main(args.dst)#fn=args.fnLOG)
    elif args.i==1:
        compare(fnLOG=args.fnLOG,fnSA=args.fnSA,fnFE=args.fnFE,ntens=args.ntens)
