id2name={'1       X' : 'H       X',
         '2       X' : 'H       X',
         '3       X' : 'H       X',
         '4       X' : 'H       X',
         '5       X' : 'H       X',
         '6       X' : 'H       X',
         '7       X' : 'H       X'}
buf=''
for line in open('peblock.pdb'):
    if 'ATOM' in line:
        key=line[13:22]
        line=line.replace(key,id2name[key])
    buf+=line  
open('peblock_allH.pdb','w').write(buf)
