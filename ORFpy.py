# A brand new simple program to get ORFs from Genome

# Uses Simple regex patterns and re.findall


import re , operator, random
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord

# Defining patterns

_start = r'ATG'
_stop = r'(?:TAG|TGA|TAA)'
_nonstop = r'(?:[CAG][TCAG]{2}|T(?:[TC][TCAG]|[AG][TC])|TGG)'
_codon = r'(?:[TCAG]{3})'
_orf_re = re.compile('(ATG'+ _nonstop + '*)\
('+ _stop +')',flags = re.I) # Edited for finding only orfs
_lead_re = re.compile(r'ATG',flags = re.I)

#define all the extractors
def extract_orfs(seqt,FROM,TO):
    orfs = []
    lens = []
    strand = []
    pos=[]
    for seqs,stra in zip([seqt, seqt.reverse_complement()],["F","RC"]):
        start_pos = []
        for a in _lead_re.finditer(str(seqs)):
            start_pos.append(a.start())
        for i in start_pos:
            a = _orf_re.search(str(seqs),i)
            if a is not None:
                s = seqt[a.start():a.end()]
                orfs.append(s)
                lens.append(len(s))
                strand.append(stra)
                if stra=="F":
                    pos.append(FROM+i)
                elif stra=="RC":
                    pos.append(TO-i)
    return (orfs,lens,strand,pos)
