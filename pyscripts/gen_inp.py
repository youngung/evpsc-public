"""
This module is to generate <vpsc7.in> file used for delta & sigma EVPSC
from command line prompt using Python
"""
import numpy as np
import argparse
import sys
import tempfile
import subprocess
import shutil
import os
import glob


num_cond_ori="""*PRECISION SETTINGS FOR CONVERGENCE PROCEDURES (default values)
1e-3 1e-3 1e-3 1e-3  ## ERRS,ERRD,ERRM,ERRSO
99  99  20    itmax:  max # of iter, external, internal and SO loops
0  2 10 2   irsvar & xrsini,xrsfin,xrstep (dummy if irsvar=0)
1              ibcinv (0: don't use <Bc>**-1, 1: use <Bc>**-1 in SC eq)
*INPUT/OUTPUT SETTINGS FOR THE RUN (default is zero)
0              irecover:read grain states from POSTMORT.IN (1) or not (0)?
0              isave:   write grain states in POSTMORT.OUT at step 'isave'?
0              icubcomp:calculate fcc rolling components?
%i              nwrite (frequency of texture downloads)
*MODELING CONDITIONS FOR THE RUN
1              interaction (0:FC,1:affine,2:secant,3:neff=10,4:tangent,5:SO)
%i  %i  %i     iupdate: update orient, grain shape, hardening
0              nneigh (0 for no neighbors, 1 for pairs, etc.)
0              iflu (0: don't calc, 1: calc fluctuations)"""


def get_num_cond(nwrite,iupori,iupsh,iuphard):
    return num_cond_ori%(nwrite,iupori,iupsh,iuphard)


## create generic dEVPSC input file
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--nph',default=1,type=int,help='# of phase')
    parser.add_argument(
        '--wph',nargs='+',default=[],help='Texture file')
    parser.add_argument(
        '--fnout',type=str,help='name of generic input file')
    parser.add_argument(
        '--eulph1',nargs='+',default=[],help='Phase 1 inclusion orientation')
    parser.add_argument(
        '--eulph2',nargs='+',default=[],help='Phase 2 inclusion orientation')
    parser.add_argument(
        '--eulph3',nargs='+',default=[],help='Phase 3 inclusion orientation')
    parser.add_argument(
        '--axiph1',nargs='+',default=[],help='Phase 1 inclusion axis')
    parser.add_argument(
        '--axiph2',nargs='+',default=[],help='Phase 2 inclusion axis')
    parser.add_argument(
        '--axiph3',nargs='+',default=[],help='Phase 3 inclusion axis')
    parser.add_argument(
        '--nwrite',default=0,type=int,help='nwrite')
    parser.add_argument(
        '--iupori',default=1,type=int,help='update orientation')
    parser.add_argument(
        '--iupsh',default=1,type=int,help='update shape')
    parser.add_argument(
        '--iuphard',default=1,type=int,help='update hardening')

    args=parser.parse_args()
    nph=args.nph
    wph=args.wph

    if len(args.eulph1)==0:args.eulph1=[0.,0.,0.]
    if len(args.eulph2)==0:args.eulph2=[0.,0.,0.]
    if len(args.eulph3)==0:args.eulph3=[0.,0.,0.]
    if len(args.axiph1)==0:args.axiph1=[1.,1.,1.]
    if len(args.axiph2)==0:args.axiph2=[1.,1.,1.]
    if len(args.axiph3)==0:args.axiph3=[1.,1.,1.]
    eulph=np.array([args.eulph1,args.eulph2,args.eulph3])
    axiph=np.array([args.axiph1,args.axiph2,args.axiph3])


    cnt=f"""1\n{nph}\t\t\tNumber of phases\n"""
    for iph in range(nph):
        cnt=f"{cnt} {wph[iph]}"
    cnt=f"{cnt}\n"
    ## phase specific inputs
    for iph in range(nph):
        cnt=f"""{cnt}** information about phase # {iph+1}\n"""
        #cnt=f"""{cnt}0\t0\t25 \n1.\t1.\t1. \n0.\t0.\t0.\n"""
        cnt=f"""{cnt}0\t0\t25\n{axiph[iph,0]}\t{axiph[iph,1]}\t{axiph[iph,2]}\n"""
        cnt=f"""{cnt}{eulph[iph,0]}\t{eulph[iph,1]}\t{eulph[iph,2]}\n"""
        cnt=f"""{cnt}*fntx\ntexture_file_generic_path_{iph}\n"""
        cnt=f"""{cnt}*fnsx\nsx_file_generic_path_{iph}\n"""
        cnt=f"""{cnt}*shape\ndummy\n"""
        cnt=f"""{cnt}*idiff\n0\ndiffraction_file_generic_path_{iph}\n"""

    num_cond=get_num_cond(args.nwrite,args.iupori,args.iupsh,args.iuphard)
    cnt=f"{cnt}{num_cond}\n"
    with open(args.fnout,'w') as fo:
        fo.write(cnt)

    if True:
        print('===')
        print(cnt)
        print('===')
