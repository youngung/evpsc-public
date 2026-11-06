import shutil
import os,tempfile,subprocess
if os.name=='nt':
    ABQ="C:/SIMULIA/EstProducts/2020/win_b64/code/bin/ABQLauncher.exe"
else:
    ABQ="/apps/commercial/abaqus/Commands/abq2020hf4"

def main(pathwork,noel,npt,instance_name,uvarmx=30,iopt=0,
         ihardlaw=0):
    """
    Arguments
    ---------
    pathwork
    noel
    npt
    instance_name
    uvarmx
    iopt (0: default; 1: vumat case)
    ihardlaw=0
    """
    import sys
    print('\n\n** Entering lookup.main')

    fn_rst=os.path.join(f"{pathwork}","tmp_%4.4i_%4.4i.rst"%(noel,npt))
    fn_odb=os.path.join(pathwork,'tmp.odb')#"%s/tmp.odb"%pathwork
    print('** <fn_rst>:',fn_rst)
    print('** <fn_odb>:',fn_odb)
    cmd = f"{ABQ} python -u abq_rst.py {noel} {npt} {fn_odb}"
    cmd = f"{cmd} {fn_rst} {uvarmx} {instance_name}"
    print('** <cmd>:',r'%s'%cmd)
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE)
    (output,err) = p.communicate()
    p_status=p.wait()
    output=output.decode('utf-8')

    with open('tmp.txt', 'w') as fo:
        fo.write(output)

    ntens=output.split("ntens:'")[1]
    ntens=int(ntens.split(')')[0].split(',')[1])
    if ntens==6: ## solid
        print('** Element type: solid')
    elif ntens==4: ## possibly shell?
        print('** Element type: possibly shell')
    else:
        msg=f'** Unexpected ntens <ntens> was given.'
        sys.stderr.write(f'{msg}\n')
        raise IOError(msg)

    if iopt==0: # SA runs only for iopt=0 case
        fn_vpsc7in="%s/vpsc7.in"%pathwork
        home=os.getcwd()
        wd=pathwork[::]
        for i in range(3):
            wd=os.path.split(wd)[0]
        print('** wd: %s'%wd)

        if os.name=='nt' and os.path.isfile(os.path.join(wd,'evpsc.exe')):
            print('** use the existing executable')
        elif os.name!='nt' and os.path.isfile(os.path.join(wd,'evpsc')):
            print('** use the existing executable')
        else:
            if os.name=='nt':
                cmd='make -j evpsc -f win_ifort.mk'
            else:
                cmd='make -j evpsc -f nix_ifort.mk'
            p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,cwd=wd)
            (output,err) = p.communicate()

        if os.name=='nt':
            fn_exe=os.path.join(wd,'evpsc.exe')
        else:
            fn_exe=os.path.join(wd,'evpsc')
        dst=pathwork[::]
        print('fn_exe:',fn_exe)

    varhist=tempfile.mktemp(prefix="youngung_tmp_")
    print('** varhist',varhist)
    cmd=f'python -u gen_hist.py -i {fn_rst} -o {varhist} -ntens {ntens}'
    print('** cmd:',cmd)
    p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
    (output,err) = p.communicate()
    p_status=p.wait()

    if iopt==0:
        print('fn_vpsc7in:',fn_vpsc7in)

        with open(fn_vpsc7in,'r') as fo:
            cnt=fo.read().split('\n')
            lines=cnt[:31]

        cnt=''
        for line in lines:
            if len(line)>0: cnt=f'{cnt}{line}\n'

        cnt=f'{cnt}*nproc \n1 \n*ivgvar(arguments) \n1 \n{varhist}'

        print('\n\n\n\n----------------')
        print('vpsc input file:')
        print(cnt)
        print('----------------\n\n\n\n')

        with open(os.path.join(dst,'vpsc7.in'),'w') as fo:
            fo.write(cnt)

        print('Path changed to dst (%s) to run evpsc-sa'%dst)
        p=subprocess.Popen(fn_exe,stdout=subprocess.PIPE,cwd=pathwork)
        (output,err)=p.communicate()
        p.wait()
        print(output.decode('utf-8'))

    cmd=f'python -u plot_abq.py -i 0 --ntens {ntens} --fnRST {fn_rst} --ihardlaw {ihardlaw}'
    if iopt==0: ## SA runs only for iopt==0
        cmd=f'{cmd} --fnSA {os.path.join(dst,"STR_STR.OUT")}'
    elif iopt==1:
        pass
    else:
        raise IOError('** Unexpected iopt given')
    if ihardlaw==20:
        cmd=f'{cmd} --fndd {os.path.join(pathwork,"DIS_DEN.OUT")}'

    print('cmd:',cmd)
    p=subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
    (output,err)=p.communicate()
    p.wait()
    print(output.decode('utf-8'))


if __name__=='__main__':
    print('** Beginning of <lookup.py>')

    import argparse, os
    import sys
    try:
        with open('run.log', 'r') as fo:
            recent = fo.read()
            recent=recent.split('\n')[-2].strip()
    except:
        print('** Could not find the most recent path from "run.log"')
        recent=None

    print('** Argument parsing')
    print('** Recent:',recent)

    with open('src/io_head.f','r') as fo:
        cnt=fo.read().lower()
        dum=cnt.split('parameter(nuvarmx=')[1].split('\n')[0]
        nuvarm=dum.split(')')[0]
        nuvarm=int(nuvarm)
        print('** <UVARM> estimated from src/io_head.f:',nuvarm)

    ## Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-p',type=str,help='path work',default=recent)
    parser.add_argument('-noel',type=int,help='abq_rst file name',default=1)
    parser.add_argument('-npt',type=int,help='dEVPSC stand-alone result',default=1)
    parser.add_argument('-instance_name',type=str,help='# of instance',default='PART-1-1')
    parser.add_argument('-nuvar',type=int,help='# of USER field var',default=nuvarm)
    parser.add_argument('-iopt',type=int,help='iopt',default=0)
    parser.add_argument('-ihardlaw',type=int,help='ihardlaw',default=0)

    args=parser.parse_args()
    main(args.p,args.noel,args.npt,args.instance_name,args.nuvar,args.iopt,args.ihardlaw)
