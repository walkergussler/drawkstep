from __future__ import absolute_import, division
import math, sys, os, time, collections, argparse, subprocess
from tempfile import NamedTemporaryFile
from itertools import izip, combinations
from Bio import SeqIO
import networkx as nx
import numpy as np

#program breaking conditions:
    #assumes same sequence will not appear twice in one file (not sure which error this will cause (maybe none?))
    #requires frequency for each read to follow last underscore in the sequence ID (continued)
    #If anything following last underscore is not frequency, you will get unexpected results
#note for fileFreqs: if f1 has 10 sequences, f2 has 20 sequences, and all of f1's seqs are in f2, then we would have fileFreqs=[10,0,10]
  
def pre_align(files): #align files using mafft before processing
    seqNum=[]
    c=0
    for file in files: 
        with open(file) as f:
            for record in SeqIO.parse(f,'fasta'):
                c+=1
        seqNum.append(c)    
    with NamedTemporaryFile() as tmpcat:
        catname=tmpcat.name
        catlist=['cat']
        for file in files:
            catlist.append(file)
        try:
            subprocess.check_call(catlist,stdout=tmpcat)
        except:
            print('Something weird happened with one of these files')
            for file in files:
                print(file)
            sys.exit()
    
        with NamedTemporaryFile() as aligned:
            alignname=aligned.name
            subprocess.check_call(['mafft', '--quiet', '--auto', '--thread', '20', '--preservecase', catname], stdout=aligned) 
            seqs=SeqIO.parse(alignname,'fasta')
            seqses=[]
            c=0
            tmp=[]
            with open(alignname) as f:
                for record in SeqIO.parse(alignname,'fasta'):
                    tmp.append(record)
                    c+=1
                    if c in seqNum:
                        seqses.append(tmp)
                        tmp=[]
    aliNames=[]
    for id in range(len(seqses)):
        seqs=seqses[id]
        print(len(seqs))
        file=files[id]
        with NamedTemporaryFile(delete=False) as f:
            SeqIO.write(seqs,f,'fasta')
            aliNames.append(f.name)
    return aliNames
    
def get_seqs(files,freqCut,output): #parse input
    print("getting seqs")
    numFiles=len(files)
    fileFreqs=np.zeros(numFiles+1)
    seqs=collections.OrderedDict()
    counter=0
    for fileNum,file in enumerate(files,start=1):
        with open(file) as input_handle:
            for record in SeqIO.parse(input_handle,"fasta"):
                if len(record.seq) > 0 and len(record.seq) < 50000:
                    try:
                        freq=int(record.id.split('_')[-1])
                    except:
                        sys.exit('something has gone wrong')
                    if record.seq in seqs: #not a new sequence
                        fileFreqs[0]+=1
                        oldsource=int(seqs[record.seq][0])
                        oldID=seqs[record.seq][1]
                        fileFreqs[oldsource]-=1
                        newfreq=int(seqs[record.seq][2]) + freq
                        seqs[record.seq]=np.array([0,oldID,newfreq])
                    else: #new sequence
                        counter+=1
                        fileFreqs[fileNum]+=1
                        seqs[record.seq]=np.array([fileNum,counter,freq])

    if freqCut>0:
        (seqs,fileFreqs)=trim_by_freq(seqs,freqCut,fileFreqs)
    print("There are %i sequences in this network" % len(seqs))
    return (seqs,fileFreqs,counter)

def trim_by_freq(seqs,freqCut,fileFreqs): #only consider sequences with frequency greater than or equal to freqCut
    print("trimming")
    newseqs=collections.OrderedDict()
    counter=0
    for seq in seqs:
        if int(seqs[seq][2])>=freqCut:
            counter+=1
            # seqs[seq][1]=str(counter)
            newseqs[seq]=[seqs[seq][0],str(counter),seqs[seq][2]]
        else:
            fileFreqs[int(seqs[seq][0])]-=1
    return (newseqs,fileFreqs)

def calc_ordered_frequencies(haploNum,haploSize,seqDict,byFreq): 
    seqs=seqDict.keys()
    freqCount = np.zeros((haploSize, 5))
    productVector = np.zeros((haploSize, 5))
    order={'A':0,'C':1,'G':2,'T':3,'-':4}
    try:
        total_reads=0
        for read in seqs:
            if byFreq:
                freq=seqs[read]
            else:
                freq=1
            total_reads+=freq
            for pos in range(haploSize):
                num=order[read[pos]]
                freqCount[pos, num] = freqCount[pos, num] + freq
        freqRel = np.divide(freqCount, float(total_reads), dtype = float)
    
    except IndexError:
        print("Your files are not aligned and it caused an error! Try again with -a")
    
    for pos in range(haploSize):
        for i in range(5):
            freqPos = freqRel[pos, i]
            if freqPos > 0:
                logFreqRel = math.log(freqPos, 2)
                productVector[pos, i] = -1*(np.multiply(freqPos, logFreqRel, dtype = float))                
    return np.sum(productVector, axis = 1)

def order_positions(hVector,seqs,haploSize): #order positions for faster building of k-step network
    invH =  np.multiply(-1, hVector, dtype = float)
    ordH = np.argsort(invH)
    # reorder the sequences by entropy
    ordSeqs = {}

    for seq in seqs:
        newOne = ''
        for p in range(haploSize):
            newOne = ''.join([newOne, seq[ordH[p]]])
        ordSeqs[newOne]=seqs.pop(seq)
        # ordSeqs.append(newOne)
    return ordSeqs

def calc_kstep(haploNum,haploSize,seqs): 
    bDict={}
    g=nx.Graph()
    t = 0
    for seq in seqs:
        [_,counter,_]=seqs[seq]
        g.add_node(counter)
        bDict[counter]=seq
    while not nx.is_connected(g):
        t = t + 1
        for seq1,seq2 in combinations(seqs.keys(),2):
            [color1,counter1,freq1]=seqs[seq1]
            [color2,counter2,freq2]=seqs[seq2]
            
            dist = 0
            for a, b in izip(seq1, seq2):
                if a != b:
                    dist+=1
                    if dist > t:
                        break
            if dist == t: 
                g.add_edge(counter1,counter2,len=dist)
    return g,bDict

def get_intermediate_sequences(seqs,g,id,bDict,haploSize): #add edges in here to be faster
    for edge1,edge2,d in g.edges_iter(data=True):
        dist=d['len']
        s1=bDict[edge1]
        s2=bDict[edge2]
        # print(dist)
        if dist!=1:
            while dist>1:
                id+=1
                for nucl in range(haploSize):
                    if s1[nucl]!=s2[nucl]:
                        takeseq1=s1[:nucl+1]+s2[nucl+1:]
                        takeseq2=s1[:nucl]+s2[nucl:]
                        if takeseq1!=s1:
                            newseq=takeseq1
                        else:
                            newseq=takeseq2
                dist-=1
                s1=newseq
                if newseq not in seqs:
                    seqs[newseq]=[8,id,1]
    return seqs

def seqs_to_nx(seqs): 
    g=nx.Graph()
    numNodes=len(seqs)
    gstr='gray{}'.format(int(math.ceil(70*(-math.exp(-numNodes/1000)+1)+20)))
        
    w_i = .07
    if numNodes <= 3000:
        w_i = math.exp(-numNodes/500)*.2+.07
    
    for seq in seqs:
        [color,id,freq]=seqs[seq]
        mult = 4
        if freq <= 2000:
            mult = -math.exp(-freq/300)*2.2+4
        w=mult*w_i
        w=str(w)
        # print(w,color,freq,w_i,mult)
        g.add_node(id,color=int(color)+1,shape='point',width=w,colorscheme='set19')
        for seq2 in seqs:
            dist=0
            
            for a, b in izip(seq, seq2):
                if a != b:
                    dist += 1
                    if dist > 1:
                        break
            if dist==1:
                g.add_edge(seqs[seq2][1],id,color=gstr)
    return g
    
MAX_LEGEND_FONT_SIZE = 200
LEGEND_FONT_SIZE = {
  100 : 6,
  500 : 10,
  1000 : 20,
  2000 : 40,
  4000 : 80,
  8000 : 160
}

def determine_legend_font_size(g): #determine font size for legend
    min_x,min_y,max_x,max_y = float('inf'),float('inf'),float('-inf'),float('-inf')
    try:
        for n,d in g.nodes_iter(data=True):
            x,y = map(float, d['pos'].split(','))
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            max_x = max(max_x, x)
            max_y = max(max_y, y)

    except KeyError:
        raise ValueError('nodes are missing the pos key, please layout graph first.')

    max_axis_len = max(max_y - min_y, max_x - min_x)
    # print('max_axis_len', max_axis_len)
    fontSize=MAX_LEGEND_FONT_SIZE
    for len_bin,font_size in sorted(list(LEGEND_FONT_SIZE.items()),key=lambda t: t[0]):
        if max_axis_len<len_bin:
            fontSize=font_size
            break
    # print('FONTSIZE:'+str(fontSize))
    return str(fontSize)

def add_legend(file,output,files,fontSize,fileFreqs,scheme): #add legend to figure
    print("Adding legend...")
    with open(file) as f:
        lines=f.readlines()
    l=len(lines)        
    lineco=0
    with open(output,"w") as g:
        for line in lines:
            lineco+=1
            if lineco==2:
                g.write('\tgraph [outputorder=edgesfirst];\n')
            if lineco!=l:
                g.write(line)
            else:      
                g.write('{  rank=sink;\n')
                g.write('\tLegend [shape=none,margin=0,pos="0,0",colorscheme='+scheme+' label=<\n')
                g.write('\t<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">\n')
                g.write('\t  <TR>\n')
                g.write('\t    <TD><FONT COLOR="%i" POINT-SIZE="%s">Shared Sequences</FONT></TD>\n' % (1,fontSize))
                g.write('\t    <TD><FONT COLOR="%i" POINT-SIZE="%s">%i</FONT></TD>\n' % (1,fontSize,fileFreqs[0]))
                g.write('\t  </TR>\n')
                
                for id in range(len(files)):
                    file=files[id]
                    shortfile=trimfh(file)
                    freq=fileFreqs[id+1]
                    g.write('\t  <TR>\n')
                    g.write('\t    <TD><FONT COLOR="%i" POINT-SIZE="%s">%s</FONT></TD>\n' % (id+2,fontSize,shortfile))
                    g.write('\t    <TD><FONT  COLOR="%i"  POINT-SIZE="%s">%i</FONT></TD>\n' % (id+2,fontSize,freq))
                    g.write('\t  </TR>\n')
                g.write('</TABLE>\n >];\n }\n}')

def trimfh(handle):
    return os.path.splitext(os.path.basename(handle))[0]

def main(prefiles,ali,output,freqCut,colorscheme,drawmode,start,keepfiles,ghost,show):
    import time
    start=time.time()
    if ali==True:
        files=pre_align(prefiles)
        file_names=prefiles
    else:
        files=prefiles
        file_names=False
    print('hello',time.time()-start)
    seqs,fileFreqs,counter=get_seqs(files,freqCut,output)
    haploNum=len(seqs)
    haploSize=len(seqs.keys()[0])
    print(len(seqs))
    # print(seqs)
    print('got seqs',time.time()-start)
    hVector=calc_ordered_frequencies(haploNum,haploSize,seqs,False)
    ordSeqs=order_positions(hVector,seqs,haploSize)
    print('prepared kstep',time.time()-start)
    print(len(ordSeqs))
    # print(ordSeqs)
    g,bDict=calc_kstep(haploNum,haploSize,ordSeqs)
    print('built kstep',time.time()-start)
    finalSeqs=get_intermediate_sequences(ordSeqs,g,counter,bDict,haploSize)
    print(len(finalSeqs))
    print('got intermediates',time.time()-start)
    # for item in finalSeqs:
        # print(finalSeqs[item])
    # g=kstep(allDistTuple,g)
    g=seqs_to_nx(finalSeqs)
    print('into nx',time.time()-start)
    tmp=nx.drawing.nx_agraph.to_agraph(g)
    if keepfiles:
        tmp.write("tmp0") 
        subprocess.check_call(['neato','-Tdot','tmp0','-o','tmp1'])
        h=nx.drawing.nx_agraph.read_dot('tmp1')
        fontSize=determine_legend_font_size(h)
        if type(file_names)==bool:
            add_legend('tmp1',output,files,fontSize,fileFreqs,colorscheme)
        else:
            add_legend('tmp1',output,file_names,fontSize,fileFreqs,colorscheme)
        print(output)
        subprocess.check_call(['neato','-Tpng',output,'-n2','-O'])
    else:
        with NamedTemporaryFile(delete=False) as tmp1:
            with NamedTemporaryFile(delete=False) as tmp0:
                tmp0name=tmp0.name
                tmp1name=tmp1.name
                tmp.write(tmp0name)
                subprocess.check_call(['neato','-Tdot',tmp0name,'-o',tmp1name])
                h=nx.drawing.nx_agraph.read_dot(tmp1name)
                fontSize=determine_legend_font_size(h)
                if type(file_names)==bool:
                    add_legend(tmp1name,output,files,fontSize,fileFreqs,colorscheme)
                else:
                    add_legend(tmp1name,output,file_names,fontSize,fileFreqs,colorscheme)
                subprocess.check_call(['neato','-Tpng',output,'-n2','-O'])
                os.unlink(output)
    print("All finished! It took",  time.time()-start)
    if show:
        try:
            subprocess.check_call(['firefox',output+'.png'])
        except:
            sys.exit()
    sys.exit()

if __name__=="__main__":
    start=time.time()
    parser=argparse.ArgumentParser(description='Show k-step network for a given pair of files')
    parser.add_argument('files',
        nargs='+',
        help="List of files to be analyzed,order of files does not matter")
    parser.add_argument('-o','--output',
        required=False,default="kstep",
        help="Your  preferred  output  name  (no  file  extension  necessary)")
    parser.add_argument('-f','--frequencycutoff',
        type=int,required=False,default=0,
        help="Minimum  frequency  sequence  must  have  to  appear  on  the  figure,will  make  the  program  run  very  slow  if passed  with -a  with >500  sequences,not reccomened  without  -a")
    parser.add_argument('-a','--align',
        action='store_true',default=False,
        help="Pass this as an argument to align your files; this will cause the program to run much more slowly if you add -f")
    parser.add_argument('-s','--show',
        action='store_true',default=True,
        help="Show picture when process is finished")
    parser.add_argument('-g','--ghost',
        action='store_true',default=False,
        help="Pass this as an argument to utilize the GHOST codebase, which could speed up performance")
    parser.add_argument('-k','--keepfiles',
        action='store_true',default=False,
        help="Pass this as an argument if you wish to keep your intermediate dot file")
    parser.add_argument('-c','--colorscheme',
        dest='graphviz_colorscheme',default='set19',
        help='The graphviz colorscheme used. There should be at least one more color than files. [default: %(default)s]')
    parser.add_argument('-d','--drawmode',
        dest='drawmode',default='allred',choices=['dont','red','white','allred','weights'],
        help="""Decide what to do with longer edges on the network. [default: %(default)s]
        #dont - remove edges from network,draw disconnected components if necessary - the disconnected components will not be placed in any specific way relative to one another by graphviz
        #allred - color longer edges red with their length set to 9 for the graphing software
        #red - only act on longer edges if they are between two different patients
        #white - all edges are drawn with their true weights,with longer links turned white
        #weights - draw all edges with labels (very busy chart),shorten longer edges""")

    args=parser.parse_args()
    files=args.files
    ali=args.align
    output=args.output
    freqCut=args.frequencycutoff
    drawmode=args.drawmode
    colorscheme=args.graphviz_colorscheme
    keepfiles=args.keepfiles
    ghost=args.ghost
    show=args.show

    if output+'.fas' in files:
        files.remove(output+'.fas')
    if len(files) > 8 and colorscheme == 'set19':
        sys.exit('The default colorscheme can only support 7 files. Please try again with a different colorscheme from this list:\nhttp://www.graphviz.org/doc/info/colors.html#brewer\nSuggestion: or just try set312')
    main(files,ali,output,freqCut,colorscheme,drawmode,start,keepfiles,ghost,show)

 
