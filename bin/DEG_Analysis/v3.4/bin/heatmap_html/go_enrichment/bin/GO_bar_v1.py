'''
This script was used to draw the most enrichment GO stat  and add hyperlinks in picture.

Author:
       huangls
Version:
        1.0;2014-12-2
'''





import sys, os, ConfigParser, argparse, glob, os.path
#from pylab import *
import numpy as np
import xlrd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import re
from operator import itemgetter, attrgetter
#sys.path.append('/share/bioCloud/huangls/02.package/python/svgwrite-1.1.6')
#import svgwrite
mpl.rcParams['interactive']=False
#from svgwrite import cm, mm 
#from cloud import join
#from casuarius import required
#from enstaller.config import default
parser = argparse.ArgumentParser(description='This script was used to draw the most enrichment GO stat  and add hyperlinks in svg picture.')
parser.add_argument('-i', '--inputFile', dest='inputFile', required=True, help='input  stat file')
parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')

parser.add_argument('-P','--pCol',dest='pCol',required=False,type=int,default=6,help='specify the P-value column,default is 5')
parser.add_argument('-O','--ontologyCol',dest='ontologyCol',required=False,type=int,default=3,help='specify the ontologyCol,default is 3')
parser.add_argument('-G','--GO_termCol',dest='GO_termCol',required=False,type=int,default=2,help='specify the GO_term column,default is 2')
parser.add_argument('-g','--geneNumberCol',dest='geneNumberCol',required=False,type=int,default=4,help='specify the geneNumber column,default 4')
parser.add_argument('-t','--topLine',dest='topLine',required=False,type=int,default=20,help='specify the topLine to draw picture,default is 20')
parser.add_argument('-X','--dataPath',dest='dataPath',required=False,default='Id path not given',help='specify the GO  id file')
parser.add_argument('-p','--outFilePrefix',dest='outFilePrefix',required=False,default='Go',help='specify the output file prefix,default is Go')
parser.add_argument('-W','--width',required=False,default=12,type=int,help='specify the width column index,default is 12')
parser.add_argument('-H','--height',required=False,default=8,type=int,help='specify the height column index,default is 8')

args = parser.parse_args()
if(not os.path.exists(args.outDir)):
    os.mkdir(args.outDir)
try:    
    assert args.inputFile
except AssertionError:
    print 'Please specify the input file "--inputFile"'
    os._exit(0)
if os.path.isfile(args.dataPath):
    args.dataPath=os.path.abspath(args.dataPath)
    
    
    
def unique(L):
    from sets import Set
    return list(Set(L))
        
        
def subStr(l,num):
    res=[]
    for s in l:
        if(len(s)>num):
            res.append(s[0:num]+' ...')
        else:
            res.append(s)
    return res
            
            
            
cmaps = [('Sequential', ['Blues', 'BuGn', 'BuPu',
'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
('Sequential (2)', ['afmhot', 'autumn', 'bone', 'cool', 'copper',
'gist_heat', 'gray', 'hot', 'pink',
'spring', 'summer', 'winter']),
('Diverging', ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
'seismic']),
('Qualitative', ['Accent', 'Dark2', 'Paired', 'Pastel1',
'Pastel2', 'Set1', 'Set2', 'Set3']),
('Miscellaneous', ['gist_earth', 'terrain', 'ocean', 'gist_stern',
'brg', 'CMRmap', 'cubehelix',
'gnuplot', 'gnuplot2', 'gist_ncar',
'nipy_spectral', 'jet', 'rainbow',
'gist_rainbow', 'hsv', 'flag', 'prism'])]




file=open(args.inputFile,'r')




data=[]
pcol=args.pCol-1
ontologyCol=args.ontologyCol-1
GO_termCol=args.GO_termCol-1
geneNumberCol=args.geneNumberCol-1
topLine=args.topLine

GO_term=[]
Ontology=[]
P_value=[]
geneNumber=[]
url=[]
for line in file:
    line.rstrip('\n')
    if line[0]=='#' or re.findall('GO_ID', line):
        continue
    tmp=re.split('\t', line)
    tmp[pcol]=float(tmp[pcol])
    data.append(tmp)
#print data
data.sort(key=lambda x:x[pcol],reverse=False)
topData=data[:topLine]
#print topData
topData.sort( key=itemgetter(ontologyCol,pcol), reverse=False)
#print topData

isff=False
for i in topData:
    f=0
    GO_term.append(i[GO_termCol])
    Ontology.append(i[ontologyCol])
    try:
        P_value.append(-math.log10(float(i[pcol])))
    except ValueError:
        P_value.append(50)
        f=1
        isff=True
    try:
        if(f==1):
            geneNumber.append(str(int(i[geneNumberCol]))+'**')
            f=0
        else:
            geneNumber.append(str(int(i[geneNumberCol])))
    except ValueError:
        if(f==1):
            tmp=int(re.split('\s+',i[geneNumberCol] )[0])
            geneNumber.append(str(tmp)+'**')
            f=0
        else:
            tmp=int(re.split('\s+',i[geneNumberCol] )[0])
            geneNumber.append(tmp)
#    url.append('javascript:window.parent.RNA_noRefReport.vennClick(\'%s\',\'%s\',\'%s\',\'%s\')' \
#               %(os.path.abspath(args.inputFile),'1',i[0],i[ontologyCol]))
    url.append('javascript:window.parent.RNA_noRefReport.DataMining.vennClick(\'%s\',\'%s\',\'%s\')' %(args.dataPath,i[0],'go'))
    
plt.style.use('ggplot')
#plt.axes(left, bottom, width,height)
GO_term=subStr(GO_term, 30)
col=[]
l=[]
cm = plt.get_cmap("Set1")
flag=0

for i,j in enumerate(Ontology):
    if (i==0):
        l.append([0,Ontology[i]])
    if(i>0 and Ontology[i] != Ontology[i-1]):
        flag+=0.125
        l.append([i,Ontology[i]])
    #print str(flag)+'\n'
    col.append(cm(flag))
    

#col=[cm(float(i)/(len(unique(Ontology)))) for i in xrange(len(unique(Ontology)))]


#fig, ax = plt.subplots(figsize=(9, 7))
#plt.subplots_adjust(left=0.115, right=0.88)

plt.figure(1, figsize=(args.width,args.height))
left,bottom,width,height=0.35, 0.1, 0.45,0.8

aa=plt.axes([left,bottom,width,height])

y_pos = np.arange(len(Ontology))
#print y_pos
#print l
mybar=aa.barh(y_pos, P_value, align='center', alpha=0.9)
for i,j in enumerate(geneNumber):
    aa.text(P_value[i]+0.01, y_pos[i], str(j))
aa.yaxis.set_ticks_position('left')
aa.yaxis.set_ticks(y_pos)
if(isff):
    xt=aa.xaxis.get_ticklocs()
    xt=[str(i) for i in xt]
    xt[-1]='oo+'
    aa.set_xticklabels(xt)
    
    
    
aa.xaxis.set_ticks_position('bottom')
aa.set_ylim([-1,len(y_pos)])
aa.set_yticklabels( GO_term)
aa.set_xlabel('-log10(p-value)')
aa.set_ylabel('Go term')
aa.set_title('The most enrichment GO terms')



for j,i in enumerate(mybar):
    
    i.set_url(url[j])
    i.set_facecolor(col[j])
bb=plt.axes([left+width+0.01, 0.1, 0.2,0.5])
try:
    bb.legend([mybar[l[2][0]],mybar[l[1][0]],mybar[l[0][0]]],[l[2][1],l[1][1],l[0][1]],loc=2)
except IndexError:
    try:
        bb.legend([mybar[l[1][0]],mybar[l[0][0]]],[l[1][1],l[0][1]],loc=2)
    except IndexError:
        bb.legend([mybar[l[0][0]]],[l[0][1]],loc=2)
        
bb.axis('off')


plt.savefig(args.outDir.rstrip('/') +'/'+args.outFilePrefix+'.svg')
plt.savefig(args.outDir.rstrip('/') +'/'+args.outFilePrefix+'.png',dpi=300)
