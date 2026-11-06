"""
Python module to archive the resulting files by changing
the extension of the resulting file names.
"""
import glob
import os
import shutil

def change(old_ext='OLD',new_ext='NEW',idel=True):
    """
    Arguments
    ---------
    old_ext
    new_ext
    idel   : Delete the old files
    """
    fns_old=glob.glob('*.%s'%old_ext)
    fns_new=[]
    for i in range(len(fns_old)):
        fn_wo_ext=fns_old[i].split('.%s'%old_ext)[0]
        fns_new.append('%s.%s'%(fn_wo_ext,new_ext))

    ## move or copy depending on <idel>
    for i in range(len(fns_old)):
        if idel:
            shutil.move(fns_old[i],fns_new[i])
            # print('%s has been replaced to %s'%(fns_old[i],fns_new[i]))
        else:
            shutil.copy(fns_old[i],fns_new[i])
            # print('%s has been copied to %s'%(fns_old[i],fns_new[i]))
    return fns_new

def main(new_ext,idel):
    fns_new1 = change('OUT',new_ext,idel)
    fns_new2 = change('out',new_ext,idel)
    fn_tot=[]
    for i in range(len(fns_new1)):
        fn_tot.append(fns_new1[i])
    for i in range(len(fns_new2)):
        fn_tot.append(fns_new2[i])
    return fn_tot

if __name__=='__main__':
    import argparse
    import os
    import sys
    import shutil
    ## Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-ext',type=str,help='New extension to use')
    parser.add_argument('-dir',type=str,default=os.path.join(os.getcwd(),'archive','tmp'),help='Directory to which the new files will be located.')
    parser.add_argument('-v',action='store_true',help="Make verbose.")
    parser.add_argument('-idel',action='store_true',help="Delete old files.")
    args        = parser.parse_args()
    ## new files
    fns = main(args.ext,args.idel)

    # print('---')
    # print('fns:')
    # for i in range(len(fns)):
    #     print(fns[i])
    # print('---')

    sys.exit()

    if len(fns)==0:
        print('No files archived.')
        os._exit(0)

    # print('**--- files created:')
    # for i in range(len(fns)):
    #     print(fns[i])

    if not(os.path.isdir(args.dir)):
        #os.mkdir(args.dir)
        os.makedirs(args.dir)

    if args.dir!=os.getcwd():
        print('**--- files copied to:')
        for i in range(len(fns)):
            # shutil.move(fns[i],args.dir)
            shutil.copy(fns[i],args.dir)
            #print(os.path.join(args.dir,fns[i]))
