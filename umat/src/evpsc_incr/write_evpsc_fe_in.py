# This module is to write evpsc_fe.in file
if __name__=='__main__':
    import argparse
    import os
    import numpy as np
    from glob import glob
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--stv',type=str,help='path to stv',
        default='c:/users/young/repo/evpsc/umat/src/evpsc_incr/temp')

    parser.add_argument(
        '--noel',nargs='+',help='(Required) Element numbers',required=True)
    parser.add_argument(
        '--npt',nargs='+',help='(Required) Integration points',required=True)

    args=parser.parse_args()

    if not os.path.exists(args.stv):
        os.makedirs(args.stv)


    ## delete stv files already there.
    fns=glob('%s/*.stv'%args.stv)
    for i in range(len(fns)):
        os.remove(fns[i])

    with open('evpsc_fe.in','w') as fo:
        fo.write('* List of elements, nodes to save (after FE simulation finishes)\n')
        fo.write('%s\n'%args.stv)

        noel=np.array(args.noel)
        npt=np.array(args.npt)
        noel=noel.astype(int)
        npt=npt.astype(int)
        # print('args.noel:',noel)
        # print('npt:',npt)
        nset=len(noel)*len(npt)

        fo.write('%i\n'%nset)
        for iel in range(len(noel)):
            for ipt in range(len(npt)):
                fo.write('%i %i\n'%(noel[iel],npt[ipt]))
