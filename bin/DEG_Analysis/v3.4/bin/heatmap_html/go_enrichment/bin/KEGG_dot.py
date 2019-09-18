'''
This script was used to draw the most enrichment KEGG stat  and add hyperlinks in picture.

Author:
       haungls
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
parser = argparse.ArgumentParser(description='This script was used to draw the most enrichment kegg stat  and add hyperlinks in picture.')
parser.add_argument('-i', '--inputFile', dest='inputFile', required=True, help='input  stat file')
parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')

parser.add_argument('-P','--pathway_termCol',dest='pathway_termCol',required=False,type=int,default=1,help='specify the pathway_term column,default is 1')
parser.add_argument('-R','--rich_factorCol',dest='rich_factorCol',required=False,type=int,default=3,help='specify the rich_factor,default is 3')
parser.add_argument('-Q','--qvalueCol',dest='qvalueCol',required=False,type=int,default=4,help='specify the qvalueCol column,default is 4')
parser.add_argument('-g','--geneNumberCol',dest='geneNumberCol',required=False,type=int,default=2,help='specify the geneNumber column,default 2')
parser.add_argument('-t','--topLine',dest='topLine',required=False,type=int,default=20,help='specify the topLine to draw picture,default is 20')
parser.add_argument('-p','--outFilePrefix',dest='outFilePrefix',required=False,default='Kegg',help='specify the output file prefix,default is Kegg')
parser.add_argument('-W','--width',required=False,default=12,type=int,help='specify the width column index,default is 12')
parser.add_argument('-X','--dataPath',required=False,default='Id path not given',help='specify the enrich gene id file')
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
pathway_termCol=args.pathway_termCol-1
rich_factorCol=args.rich_factorCol-1
qvalueCol=args.qvalueCol-1
geneNumberCol=args.geneNumberCol-1
topLine=args.topLine

pathway_term=[]
rich_factor=[]
qvalue=[]
geneNumber=[]
geneID=[]
url=[]
for line in file:
    line.rstrip('\n')
    
    if line[0]=='#' or re.findall('pathway', line):
        continue
    tmp=re.split('\t', line)
    tmp[rich_factorCol]=float(tmp[rich_factorCol])
    tmp[qvalueCol]=float(tmp[qvalueCol])
    
    data.append(tmp)
#print data
data.sort(key=lambda x:x[qvalueCol],reverse=False)
topData=data[:topLine]

topData.sort( key=itemgetter(qvalueCol,pathway_termCol), reverse=False)
#print topData
for i in topData:
    pathway_term.append(i[pathway_termCol])
    rich_factor.append(float(i[rich_factorCol]))
    qvalue.append(-math.log10(float(i[qvalueCol])))
    geneID.append(i[5])
    try:
        geneNumber.append(int(i[geneNumberCol]))
    except ValueError:
        tmp=int(re.split('\s+',i[geneNumberCol] )[0])
        geneNumber.append(tmp)
    #for j in i:
    #    if ";" in j:
    #        url.append('javascript:window.parent.RNA_noRefReport.vennClick(\'%s\',\'%s\',\'%s\')'%(args.dataPath,i[7],'kegg'))
    
    try:
        url.append('javascript:window.parent.RNA_noRefReport.DataMining.vennClick(\'%s\',\'%s\',\'%s\')'\
                %(args.dataPath,i[7],'kegg'))
    except IndexError:
        url.append('javascript:window.parent.RNA_noRefReport.DataMining.vennClick(\'%s\',\'%s\',\'%s\')'\
                %(args.dataPath,i[1],'kegg'))

#fig, ax = plt.subplots(figsize=(9, 7))
#plt.subplots_adjust(left=0.115, right=0.88)

#plt.style.use('ggplot')
#plt.axes(left, bottom, width,height)





left,bottom,width,height=0.35,0.1,0.5,0.8
plt.figure(1, figsize=(args.width,args.height))
aa=plt.axes([left,bottom,width,height])
pathway_term=subStr(pathway_term, 55)
y_pos = np.arange(len(pathway_term))
#print y_pos
#print qvalue
#mybar=aa.barh(y_pos, P_value, align='center', alpha=0.9)
dot=10
ec=aa.scatter(rich_factor, y_pos, s=np.array(geneNumber)*dot, c=qvalue, alpha=1)
ec.set_urls(url)
#aa.yaxis.set_ticks_position('left')
aa.yaxis.set_ticks(y_pos)
#aa.xaxis.set_ticks_position('bottom')
aa.set_ylim([-1,len(y_pos)])
aa.set_yticklabels( pathway_term)
aa.set_xlabel('Rich factor')
#aa.set_ylabel('Go term')
aa.set_title('Statistics of Pathway Enrichment')
aa.grid(True)

#ec.set_array(np.array(qvalue).ravel())
#cbar = plt.colorbar(ec)
#cbar.set_label('q-value')
#print cbar

bb=plt.axes([left+width+width/11, height*0.6 +bottom, width*0.04,height*0.3])
cmap = mpl.cm.jet

norm = mpl.colors.Normalize(vmin=min(qvalue), vmax=max(qvalue))
cb1 = mpl.colorbar.ColorbarBase(bb, cmap=cmap,
                                   norm=norm,
                                   orientation='vertical')
bb.set_title('q-value')

cc=plt.axes([left+width+width/12, height*0.15 +bottom, width*0.1,height*0.35])

dotsize=np.arange(math.ceil(max(geneNumber)/10)*10,0,-10)

cc.scatter([1]*len(dotsize),np.arange(1,len(dotsize)+1,1),c='k',s=dotsize*dot)
cc.set_xlim(0,2)
cc.set_ylim(0,len(dotsize)+0.5)
for i,j in enumerate(dotsize):
    cc.text(1.7,i+0.9,str(int(j)))
cc.set_title('Gene Number')
cc.axis('off')
                                    
# for j,i in enumerate(mybar):
#     
#     i.set_url(url[j])
#     i.set_facecolor(col[j])
# bb=plt.axes([0.8, 0.1, 0.2,0.5])
# bb.legend([mybar[l[2][0]],mybar[l[1][0]],mybar[l[0][0]]],[l[2][1],l[1][1],l[0][1]],loc=2)
# bb.axis('off')


plt.savefig(args.outDir.rstrip('/') +'/'+args.outFilePrefix+'.svg')
plt.savefig(args.outDir.rstrip('/') +'/'+args.outFilePrefix+'.png',dpi=300)



    
















