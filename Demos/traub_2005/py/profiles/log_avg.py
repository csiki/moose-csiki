import os
import re
from numpy import *

r = re.compile('[ \t\n\r]+')
r.split("abc:def  ghi")


funcdic = dict()
for fname in os.listdir("."):
    if fname.endswith('.log'):
        with open(fname) as f:
            for line in f:
                linelist = r.split(line)[1:-1]
                for i in range(5):
                    if linelist[i][-1] == '%':
                        linelist[i] = linelist[i][:-1]
                    linelist[i] = float(linelist[i])
                
                linelist.insert(0, 1.0)
                key = linelist[-1]
                linelist = linelist[:-1]
                if key in funcdic:
                    funcdic[key] = funcdic[key] + array(linelist)
                else:
                    funcdic[key] = array(linelist)

# avg
for data_arr in funcdic.values():
    data_arr /= data_arr[0]

for key in funcdic:
    for val in funcdic[key]:
        print '%0.2f' % val, '\t\t',
    print key
