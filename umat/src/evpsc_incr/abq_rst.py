"""
How to use this module:

- Run example:
c:\> abaqus python -u abq_rst.py <noel> <npt> <odb file name> <output file> <uvarm_num> <instance_name>

- Depending on the number of arguments, the behavior of this Python script for Abaqus changes.
** case 1, without any argument at all.
  c:\> abaqus python abq_rst.py
   In this case, abq_rst.py assumes using 'noel of 0', 'npt of 1',
   and 'the latest odb' from 'run.log'. The UVARM_num is set to 30, and instance_name is set to 'part-1-1'.
   The results are written to 'abq_rst.txt'
** case 2, with 2 arguments (<noel>, <npt>,)
   the given arguments in that order are used.
** case 3, with 5 arguments (<noel>, <npt>, <odb filename>, <output file>, <uvarm_num>)
   the given arguments in that order are used.
** case 4, with 6 arguments (<noel>, <npt>, <odb filename>, <output file>, <uvarm_num>, <instance_name>)
   the given arguments in that order are used.

In the case of vumat.f, UVARM is not available, so that <uvarm_num> should be 0
"""
print("** Getting started with Abaqus Python script to extract element field output (COORD, S, LE, ER, UVARM) as function of time **")
import sys
import os

print 'len(sys.argv)==',len(sys.argv)

if len(sys.argv)==1: # Case 1: no cargument at all
    with open('run.log','r') as fo:
        cnt=fo.read()
        last_loc=cnt.split('\n')[-1]
        print('lost loc:',last_loc)
        print('len(sys.argv)',len(sys.argv))
    noel=0
    npt=1
    fnin=os.path.join(last_loc,'tmp.odb')
    fnout='abq_rst.txt'
    UVARM_num=30
    instance_name='part-1-1'

elif len(sys.argv)==3: # Case 2: two rguments (<noel>, <npt>)
    noel =int(sys.argv[1])-1
    npt =int(sys.argv[2])

    with open('run.log','r') as fo:
        cnt=fo.read()
        last_loc=cnt.split('\n')[-1]
        print('lost loc:',last_loc)
        print('len(sys.argv)',len(sys.argv))
    fnin=os.path.join(last_loc,'tmp.odb')
    fnout='abq_rst.txt'
    UVARM_num=30
    instance_name='part-1-1'
elif len(sys.argv)==6: # Case 3: 5 arguments
    noel =int(sys.argv[1])-1
    npt =int(sys.argv[2])
    fnin=sys.argv[3]
    fnout=sys.argv[4]
    UVARM_num=int(sys.argv[5])
    instance_name='part-1-1'
elif len(sys.argv)==7: # Case 4: 6 arguments
    noel =int(sys.argv[1])-1
    npt =int(sys.argv[2])
    fnin=sys.argv[3]
    fnout=sys.argv[4]
    UVARM_num=int(sys.argv[5])
    instance_name=sys.argv[6]
else:
    print('** Error: the number of argument is not what is expected')
    sys.exit(-1)


import glob
if not(os.path.isfile(fnin)):
    print '** could not find %s'%fnin
    fnins=glob.glob(fnin)
    fnins.sort()
    fnin=fnins[-1]
    if not(os.path.isfile(fnin)):
        print '** there is no %s'%fnin
        sys.exit(-1)


from odbAccess import *
from abaqusConstants import *

from odbMaterial import *
from odbSection import *


odb = openOdb(path=fnin)
myAssembly = odb.rootAssembly

odb.rootAssembly.instances.keys()
print('odb.rootAssembly.instances:', odb.rootAssembly.instances)
print('odb.rootAssembly.instances.keys():', odb.rootAssembly.instances.keys())
for i, key in enumerate(odb.rootAssembly.instances.keys()):
    print(key)
print'-------------'
element_arrays=odb.rootAssembly.instances[instance_name.upper()].elements  ## element arrays

step1 = odb.steps.values()[0]
try:
    StrainField=step1.frames[0].fieldOutputs['LE']
except:
    StrainField=step1.frames[0].fieldOutputs['E']

print('----------------------------------')
print('Keys available in fieldoutputs')
for key in step1.frames[0].fieldOutputs.keys():
    print(key)

print('----------------------------------')



print('len(element_arrays):',len(element_arrays))
e1=element_arrays[noel] # 1st element
field=StrainField.getSubset(region=e1,position=INTEGRATION_POINT)

ntens=0
for v in field.values:
    if v.integrationPoint==npt:
        ntens=len(v.data)

print('** ntens:',ntens)
print('** npt  :',npt)

fo = open(fnout,'w')
head=['step','incr','time']
fmt=['%5s','%5s','%20s']
for i in range(3):
    fmt.append('%20s') ## LE
    head.append('%s'%(['X','Y','Z'][i]))
for i in range(ntens):  ## LE
    fmt.append('%20s')
    head.append('LE%i'%(i+1))
for i in range(ntens):  ## S
    fmt.append('%20s')
    head.append('S%i'%(i+1))
for i in range(ntens):  ## ER
    fmt.append('%20s')
    head.append('ER%i'%(i+1))
for i in range(UVARM_num):  ## UVARM
    fmt.append('%20s')
    head.append('UVARM%i'%(i+1))

#print fmt
#print head
for i in range(len(head)):
    fo.write(fmt[i]%head[i])
fo.write('\n')


t0=0.
kinc_t=0
for istp in range(len(odb.steps.values())):
    print('step:%i'%(istp+1))
    step1 = odb.steps.values()[istp]
    nframe=len(step1.frames)
    element_arrays=odb.rootAssembly.instances[instance_name.upper()].elements  ## element arrays
    print('nframe:',nframe)

    for incr in range(nframe):
        # print('incr:',incr)
        kinc_t=kinc_t+1
        Coords=step1.frames[incr].fieldOutputs['COORD']

        try:
            StrainField=step1.frames[incr].fieldOutputs['LE']
        except:
            StrainField=step1.frames[incr].fieldOutputs['E']
        StressField=step1.frames[incr].fieldOutputs['S']
        SRField=step1.frames[incr].fieldOutputs['ER']


        if istp==0 and incr==0:
            print('StrainField')
            print(StrainField)
            print('StressField')
            print(StressField)
            print('SRField')
            print(SRField)


        #odb.rootAssembly.instances.values()[0].elements
        t=step1.frames[incr].frameValue + t0
        if incr==nframe-1: t0=t
        e1=element_arrays[noel] # 1st element

        field=Coords.getSubset(region=e1,position=INTEGRATION_POINT) # from element array
        k=0
        # print('field.values',field.values)
        for v in field.values:
            # print('v:',v)
            if v.integrationPoint==npt: ## integration point
                # print('v.integrationPoint:',v.integrationPoint)
                # print('npt:',npt)
                fo.write('%5.4i'%istp)
                fo.write('%5.4i'%incr)
                fo.write('%20.12e'%t)
                # print('len(v.data):',len(v.data))
                for i in range(len(v.data)):
                    k=k+1
                    fo.write('%+20.12e'%v.data[i])

        # print('COORD len(v.data)):', len(v.data))
        # print('** k',k)

        #if istp==0 and incr==20: print('coords has %i components'%len(v.data))

        field=StrainField.getSubset(region=e1,position=INTEGRATION_POINT)


        for v in field.values:
            # print(' v.elementLabel:',  v.elementLabel)
            # print('v.integrationPoint:',v.integrationPoint)
            if v.integrationPoint==npt: ## integration point
                for i in range(len(v.data)):
                    k=k+1
                    # print(v.data[i])
                    fo.write('%+20.12e'%v.data[i])
                break

        # print('StrainField len(v.data)):', len(v.data))
        # print('** k',k)

        if istp==0 and incr==0: print('strain has %i components'%len(v.data))

        field=StressField.getSubset(region=e1,position=INTEGRATION_POINT)
        for v in field.values:
            # print v.elementLabel
            if v.integrationPoint==npt: ## integration point
                for i in range(len(v.data)):
                    fo.write('%+20.12e'%v.data[i])
                break
        # print('StressField len(v.data)):', len(v.data))

        if istp==0 and incr==20: print('stress has %i components'%len(v.data))

        field=SRField.getSubset(region=e1,position=INTEGRATION_POINT)
        for v in field.values:
            # print v.elementLabel
            if v.integrationPoint==npt: ## integration point
                for i in range(len(v.data)):
                    fo.write('%+20.12e'%v.data[i])
                break

        #print('SRField len(v.data)):', len(v.data))

        if istp==0 and incr==20: print('SR has %i components'%len(v.data))

#-------------------------------------------------------------------------------
        if incr==0:# and False:
            print('uvarm_num:',UVARM_num)
            print('keys:',step1.frames[incr].fieldOutputs.keys())

        #try:

        for i in range(UVARM_num):
            UVARM=step1.frames[incr].fieldOutputs['UVARM%i'%(i+1)]
            field=UVARM.getSubset(region=e1,position=INTEGRATION_POINT)
            for v in field.values:
                if v.integrationPoint==npt: ## integration point
                    fo.write('%+20.12e'%v.data)
                    break


            # if istp==0 and incr==20: print('uvar %i components'%(i+1))

        #except:
        #    print('** Warning: could not obtain uvarm...')
        #    pass

        # if istp==1:
        #     print("Stopped for debugging")
        #     sys.exit(0)

        fo.write('\n')

print('kinc_t:',kinc_t)
fo.close()
