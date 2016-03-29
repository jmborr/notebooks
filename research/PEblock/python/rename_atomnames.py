id2name={'1       X' : 'PN      X',
         '2       X' : 'PH      X',
         '3       X' : 'PC      X',
         '4       X' : 'H       X',
         '5       X' : 'T       X',
         '6       X' : 'N       X',
         '7       X' : 'P       X'}
buf=''
for line in open('peblock.pdb'):
    if 'ATOM' in line:
        key=line[13:22]
        line=line.replace(key,id2name[key])
    buf+=line  
open('peblock_amber.pdb','w').write(buf)
