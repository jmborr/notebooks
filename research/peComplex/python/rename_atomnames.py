id2name={'1       X':'PN      X',
'2       X':'PC      X',
'3       X':'H       X',
'4       X':'T       X',
'5       X':'N       X',
'6       X':'P       X',
}
buf=''
for line in open('pecomplex.pdb'):
    if 'ATOM' in line:
        key=line[13:22]
        line=line.replace(key,id2name[key])
    buf+=line  
open('pecomplex_amber.pdb','w').write(buf)
