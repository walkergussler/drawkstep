from __future__ import absolute_import
from __future__ import division
import math, sys, os, re, time, collections, argparse, copy, subprocess
from scipy.sparse.csgraph import connected_components,csgraph_from_dense
from tempfile import NamedTemporaryFile
from itertools import izip,combinations,chain
from Bio import SeqIO
import networkx as nx
import numpy as np
from ghost.util.distance import hamming

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
        with open(file,"rU") as input_handle:
            for record in SeqIO.parse(input_handle,"fasta"):
                if len(record.seq) > 0 and len(record.seq) < 50000:
                    try:
                        freq=int(record.id.split('_')[-1])
                    except:
                        print(record.id)
                        raw_input()
                    if record.seq in seqs: #not a new sequence
                        fileFreqs[0]+=1
                        oldsource=int(seqs[record.seq][0])
                        oldID=seqs[record.seq][1]
                        fileFreqs[oldsource]-=1
                        newfreq=int(seqs[record.seq][2]) + freq
                        seqs[record.seq]=np.array(map(str,['0',oldID,newfreq]))
                    else: #new sequence
                        counter+=1
                        fileFreqs[fileNum]+=1
                        seqs[record.seq]=np.array(map(str,[fileNum,counter,freq]))

    if freqCut>0:
        (seqs,fileFreqs)=trim_by_freq(seqs,freqCut,fileFreqs)
    print("There are %i sequences in this network" % len(seqs))
    return (seqs,fileFreqs)

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

def make_allDistTuple(seqs,scheme): #compute distances (manually) between all sequences, initialize graph with nodes
    print("Making kstep network...")
    ourStructure=[]

    for item in combinations(seqs,2):
        seq1=item[0]
        seq2=item[1]
        dist=sum(0 if a == b else 1 for a,b in izip(seq1,seq2))
        node1=int(seqs[seq1][1])
        node2=int(seqs[seq2][1])
        ourStructure.append((dist,(node1,node2)))
        retstruct=iter(sorted(ourStructure,key=lambda t:t[0]))

    g=nx.Graph()
    for seq in seqs:
        cstr=str(int(seqs[seq][0])+1)
        g.add_node(int(seqs[seq][1]),colorscheme=scheme,color=cstr,shape='point',freq=seqs[seq][2])
    return retstruct,g
        
def make_allDistTuple_fast(seqs,scheme): #compute distances (with ghost's hamming function) between all sequences, initialize graph with nodes
        print("Making kstep network (with ghost)")
        dist_array=calc_distance_matrix(seqs)
        ourStructure=[]
        for item in combinations(seqs,2):
                seq1=item[0]
                seq2=item[1]
                node1=int(seqs[seq1][1])
                node2=int(seqs[seq2][1])                
                dist=dist_array[node1-1,node2-1]
                ourStructure.append((dist,(node1,node2)))
        retstruct=iter(sorted(ourStructure,key=lambda t:t[0]))
        g=nx.Graph()
        for seq in seqs:
                cstr=str(int(seqs[seq][0])+1)
                g.add_node(int(seqs[seq][1]),colorscheme=scheme,color=cstr,shape='point',freq=seqs[seq][2])
        return retstruct,g
        
def calc_distance_matrix(finalSeqs): #calculate distance matrix from list of sequences using ghost's hamming function
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr

def kstep(distances,g,threshold=float('inf')):
    """\
    Build a k-step network.
    :param d_iter: An iterable which yields distances in the form of a tuple of tuples (distance,(node name 1,node name 2))
    :param g: A NetworkX graph initialized with nodes using ids from the iterable.
    :param connect: Throw an exception if the graph is not connected by the k-step network.
    """
    d_iter=iter(distances)
    uf=union_find(g.nodes())
    current,(i,j)=next(d_iter)
    d_next,(n_i,n_j)=current,(i,j)
    
    try:
        while current<threshold and not uf.connected():
            next_uf=copy.deepcopy(uf)
            
            while d_next == current:
                if uf.find(n_i) != uf.find(n_j):
                    g.add_edge(n_i,n_j,len=d_next)
                    next_uf.join(n_i,n_j)
                d_next,(n_i,n_j)=next(d_iter)

            uf=next_uf
            current=d_next
            i=n_i
            j=n_j

    except StopIteration:
        if current<threshold and not uf.connected():
            raise ValueError('kstep network did not connect network.')
    return g

class union_find(object): # this could be replaced with an import statement
    ### An implementation of union find data structure.
    ### It uses weighted quick union by rank with path compression.

    def __init__(self,node_ids):
        """\
        Initialize an empty union find object.
        :param node_ids: The identifiers of the nodes.
        """
        self._sets={
                node_id : {
                    'rank' : 0,
                    'parent' : node_id
                }  for idx,node_id in enumerate(node_ids)
        }
        self.component_count=len(self._sets)

    def find(self,x):
        try:
            p_idx=self._sets[x]['parent']
            if p_idx != x:
                self._sets[x]['parent']=self.find(self._sets[p_idx]['parent'])
            return self._sets[x]['parent']
        except KeyError:
            raise KeyError('ID {0}is not a member of the union'.format(x))

    def join(self,p,q):
        # # # Combine sets containing p and q into a single set.
        p_id=self.find(p)
        q_id=self.find(q)
        pRoot=self._sets[p_id]
        qRoot=self._sets[q_id]
        if p_id != q_id:
            self.component_count -= 1
        if pRoot['rank']<qRoot['rank']:
            pRoot['parent']=q_id
        elif pRoot['rank'] > qRoot['rank']:
            qRoot['parent']=p_id
        else:
            qRoot['parent']=p_id
            pRoot['rank'] += 1

    def connected(self):
        it=iter(self._sets)
        f=self.find(next(it))
        for i in it:
            if self.find(i) != f:
                return False
        return True

def edit_graph(g,numNodes,drawmode): #dynamically determine some properties of the figure such as edge color and node size
    print("Editing graph...")
    gstr='gray{}'.format(int(math.ceil(70*(-math.exp(-numNodes/1000)+1)+20)))
    w_i = .07
    if numNodes <= 3000:
        w_i = math.exp(-numNodes/500)*.2+.07
    for item in g.nodes():
        freq=int(g.node[item]['freq'])
        mult = 4
        if freq <= 2000:
            mult = -math.exp(-freq/300)*2.2+4
        w=mult*w_i
        w=str(w)
        g.node[item]['width']=w

    for u,v,data in g.edges(data=True):
        dist=int(data['len'])
        if dist<10:
            data['color']=gstr
            if drawmode=='weights':
                data['label']=dist
                data['fontsize']=24
        else:
            if drawmode=='dont': # we probably dont need these lines
                g.remove_edge(u,v) #we probably dont need these lines
            elif drawmode=='white':
                color1=g.node[u]['color']
                color2=g.node[v]['color']
                if color1!=color2 and color1!=0 and color2 !=0:
                    data['color']='white'
                    data['len']=9
            elif drawmode=='red':
                color1=g.node[u]['color']
                color2=g.node[v]['color']
                if color1!=color2 and color1!=0 and color2 !=0:
                    data['color']='red'
                    data['len']=9
                    data['label']=dist
                    data['fontsize']=24
            elif drawmode=='allred':
                data['color']='red'
                data['len']=9
                data['label']=dist
                data['fontsize']=24
            elif drawmode=='weights':
                data['color']=gstr
                data['label']=dist
                data['len']=9
                data['fontsize']=24
            else:
                sys.exit("You have entered an invalid option for drawmode. Look at help for more info")
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

def main(prefiles,ali,output,freqCut,colorscheme,drawmode,start,keepfiles,ghost):
    if ali==True:
        files=pre_align(prefiles)
        file_names=prefiles
    else:
        files=prefiles
        file_names=False
    
    seqs,fileFreqs=get_seqs(files,freqCut,output)
    if ghost:
        allDistTuple,g=make_allDistTuple_fast(seqs,colorscheme)
    else:
        allDistTuple,g=make_allDistTuple(seqs,colorscheme)        
    
    if drawmode=='dont':
        g=kstep(allDistTuple,g,threshold=10)
    else:
        g=kstep(allDistTuple,g)

    numNodes=sum(fileFreqs)
    g=edit_graph(g,numNodes,drawmode)
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
    end=time.time()
    diff=end-start
    print("All finished! It took %.2f seconds to produce your chart" % diff)
    ##########uncomment this part if you want to display your image as soon as you're done processing
    # try:
        # subprocess.check_call(['firefox',output+'.png'])
    # except:
        # sys.exit()
    ##########
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

    if output+'.fas' in files:
        files.remove(output+'.fas')
    if len(files) > 8 and colorscheme == 'set19':
        sys.exit('The default colorscheme can only support 7 files. Please try again with a different colorscheme from this list:\nhttp://www.graphviz.org/doc/info/colors.html#brewer\nSuggestion: or just try set312')
    main(files,ali,output,freqCut,colorscheme,drawmode,start,keepfiles,ghost)

 