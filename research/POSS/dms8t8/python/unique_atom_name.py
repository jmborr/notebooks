#!/usr/env/python

from pdb import set_trace as tr
import os
import re
wd='/SNSlocal/projects/jbq/POSS/dms8t8'
inpdb='{0}/dms8t8.pdb'.format(wd)
outpdb='{0}/dms8t8_nd.pdb'.format(wd)

last_indexes={}
buf=''
for line in open(inpdb,'r').readlines():
    if line[0:5]=='ATOM ':
        name=line[12:16]
        element=re.search('\s*([a-zA-Z]+)\d+\s*',name).group(1)  #Example: 'Si10' extracts 'Si'
        last_index='1 '
        if element not in last_indexes.keys():
            last_indexes[element]=1
        else:
            last_indexes[element] += 1
            last_index='{0:<2d}'.format(last_indexes[element]) #right padding with spaces
        newname=re.sub('\d+\s*', last_index, name) #replace indexes and right padding spaces
        line = re.sub(name, newname, line)
    buf+=line

open(outpdb,'w').write(buf)
