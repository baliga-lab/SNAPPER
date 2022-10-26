#!/usr/bin/python3.7
"""@Author: Vivek Srinivas
@Affiliation: Baliga Lab, ISB, Seattle, WA, USA
@Title: SNAPPER
@Version: 1.0
"""

import os, re, pickle
from Bio import SeqIO
import pandas as pd
import ORFpy


#helper classes

class RangeDict(dict):
    def __getitem__(self, item):
        if not isinstance(item, range): 
            for key in self:
                if item in key:
                    return self[key]
            raise KeyError(item)
        else:
            return super().__getitem__(item)

# main class

class aligner:
    def __init__(self,reference_genome,annotation_file,sample_name,*reads):
        self.ref = reference_genome
        self.index_name = reference_genome.split("/")[-1].split(".")[0]
        self.reads=reads
        self.sample_check = None
        self.annotations = pd.read_csv(annotation_file,delimiter="\t",index_col=0)
        self.annotations_pickled_dict = "%s_annotations_dict.pickle"%(annotation_file.split(".")[0])
        self.features_pickled_dict = "%s_features_dict.pickle"%(annotation_file.split(".")[0])
        self.annotations_dict = None
        self.feature_dict = None
        self.syn_pos = []
        self.ns_pos = []
        self.reference_length = len(list(SeqIO.parse(reference_genome,"fasta"))[0].seq)
        self.sample_name = sample_name
        self.parent_directory = os.getcwd()
        self.out_directory = None
        self.bam_bowtie = None#
        self.bam_bwa = None#
        self.bam_combined = None#
        self.qual=20
        self.variants_bcf = None#
        self.variants_fb = None#
        self.variants_combined = None#
        self.variants_annotated = None#
        self.clean = True
        self.bowtie_coverage = None
        self.bwa_coverage = None

    
    #create output directory
    def create_out_directory(self,out_directory=None):
        if out_directory is None:
            self.out_directory = self.parent_directory+"/"+self.sample_name+"_SF"
        else:
            self.out_directory = out_directory
        if os.path.isdir(self.out_directory):
            pass
        else:
            os.mkdir(self.out_directory)


    #concatenate fasta
    def check_samples(self):
        if len(self.reads) ==2:
            pass
        elif len(self.reads) > 2 and len(self.reads)%2==0:
            r_dict = pd.DataFrame(columns=["Lane","Read","File"])
            for n,read in enumerate(self.reads):
                lane,pair = read.split("_")[-3:-1]
                r_dict.at[n,["Lane","Read","File"]]= lane,pair,read
            r_dict = r_dict.sort_values(by=['Read', 'Lane'])
            cat_reads = []
            for gn,group in r_dict.groupby("Read"):
                out_dir = "%s/%s_cat"%(self.parent_directory,self.sample_name)
                out_file = out_dir+"/"+"%s_%s.fastq"%(self.sample_name,gn)
                cat_reads.append(out_file)
                if os.path.isfile(out_file):
                    pass
                else:
                    cat_cmd = "cat %s > %s"%(" ".join(group.File),out_file)
                    if os.path.isdir(out_dir):
                        pass
                    else:
                        os.mkdir(out_dir)
                    os.system(cat_cmd)
            self.reads = cat_reads
        else:
            print("Reads data (fasta) should be in pairs")
    
    #indexing functions

    def create_index_bowtie(self):
        r = re.compile(self.index_name+".*bt2")
        indexes = list(filter(r.match, os.listdir("references_and_indexes")))
        if len(indexes) == 6:
            pass
        else:
            indexer_command = "bowtie2-build %s %s"%(self.ref,self.index_name)
            os.system(indexer_command)
            print("bowtie2 indexes are created in 'references_and_indexes' directory")

    def create_index_bwa(self):
        r = re.compile(self.index_name+".*(?:amb|ann|bwt|pac|sa)")
        indexes = list(filter(r.match, os.listdir("references_and_indexes")))
        if len(indexes) == 5:
            pass
        else:
            indexer_command = "bwa index %s"%(self.ref)
            os.system(indexer_command)
            print("bwa indexes are created in 'references_and_indexes' directory")

    #function to check coverage

    def check_coverage(self, alignment):
        if os.path.exists("%s.bai"%alignment):
            pass
        else:
            os.system("samtools index %s"%alignment)
        coverage = int(os.popen('samtools mpileup %s | awk -v X="${MIN_COVERAGE_DEPTH}"'%(alignment)+" '$4>=X' | wc -l").read())
        return coverage/self.reference_length

    #functions to align genome sequence
            
    def align_reads_bowtie(self):
        _input = ",".join(self.reads)
        self.bam_bowtie = self.out_directory+"/"+self.sample_name+"_bowtie.bam"
        if os.path.isfile(self.bam_bowtie):
            pass
        else:
            bowtie_command = "bowtie2 -x %s --no-unal -U %s -S - -p 12 | samtools view -bS - | samtools sort -m 5G -o %s"%(self.ref.split(".fasta")[0],_input,self.bam_bowtie)
            #print(bowtie_command)
            os.system(bowtie_command)
            self.bowtie_coverage = self.check_coverage(self.bam_bowtie)           

    def align_reads_bwa(self):
        _input = " ".join(self.reads)
        self.bam_bwa = self.out_directory+"/"+self.sample_name+"_bwa.bam"
        if os.path.isfile(self.bam_bwa):
            pass
        else:
            bwa_command = "bwa mem %s %s | samtools sort -o %s"%(self.ref,_input,self.bam_bwa)
            #print(bwa_command)
            os.system(bwa_command)
            self.bwa_coverage = self.check_coverage(self.bam_bwa)

    # function to combine alignments and produce concensus

    def merge_alignments(self):
        self.bam_combined = self.out_directory+"/"+self.sample_name+"_combined.bam"
        if os.path.isfile(self.bam_combined):
            pass
        else:        
            merge_command = "samtools merge %s %s %s"%(self.bam_combined,self.bam_bowtie,self.bam_bwa)
            os.system(merge_command)

    # functions to call variants

    def call_variants_bcftools(self,alignment = None):
        self.variants_bcf = self.out_directory+"/"+self.sample_name+"_bcf.vcf"
        if os.path.isfile(self.variants_bcf):
            pass
        else:    
            if alignment is None:
                bcftools_command = "bcftools mpileup -Ou -f %s %s | bcftools call -mv --ploidy 1  -Ov -o %s"%(self.ref,self.bam_combined,self.variants_bcf)
            else:
                bcftools_command = "bcftools mpileup -Ou -f %s %s | bcftools call -mv --ploidy 1  -Ov -o %s"%(self.ref,alignment,self.variants_bcf)
            os.system(bcftools_command)

    def call_variants_freebayes(self,alignment = None):
        self.variants_fb = self.out_directory+"/"+self.sample_name+"_fb.vcf"
        if os.path.isfile(self.variants_fb):
            pass
        else:  
            if alignment is None:
                fb_command = 'freebayes -f %s -b %s -v %s -p 1 -q 20'%(self.ref,self.bam_combined,self.variants_fb)
            else:
                fb_command = 'freebayes -f %s -b %s -v %s -p 1 -q 20'%(self.ref,alignment,self.variants_fb)
            os.system(fb_command)

    # functions to filter and combine variants
    def f_and_c_variants(self):
        self.variants_combined = self.out_directory+"/"+self.sample_name+"_combined_variants.csv"
        if os.path.isfile(self.variants_combined):
            pass
        else:  
            comb_data = []
            for calls,caller,start in zip([self.variants_bcf, self.variants_fb],["BCF","FB"],[28,60]):
                data = pd.read_csv(calls,sep="\t",skiprows=start,index_col=1)[["QUAL","REF","ALT"]]
                data["CALLER"] = caller
                comb_data.append(data[data["QUAL"]>=self.qual])
                
            for i,rows in comb_data[1].iterrows():
                if i in comb_data[0].index:
                    comb_data[0].at[i,"CALLER"]="BCF,FB"
                else:
                    comb_data[0].loc[i]=rows

            comb_data[0].to_csv(self.variants_combined)
        
    
    # functions to annotate variants

    def chunks(self,lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def get_syn_ns(self,p1,p2):
        syn_pos, non_syn_pos = [],[]
        if p1 < p2:
            r = list(range(p1,p2))
        elif p1>p2:
            r = list(reversed(range(p1,p2)))

        syn_pos,non_syn_pos = [],[]

        for i in self.chunks(r,3):
            self.ns_pos.extend(i[:2])
            self.syn_pos.extend([i[2]])
            

    def create_annotations_dict(self):


        if os.path.isfile(self.annotations_pickled_dict) and os.path.isfile(self.features_pickled_dict):
            with open(self.annotations_pickled_dict,"rb") as handle1:
                self.annotations_dict = pickle.load(handle1)
                
            with open(self.features_pickled_dict,"rb") as handle2:
                self.feature_dict = pickle.load(handle2)

        else:
            #create dict
            gene_dict = {}
            feature_dict = {}
           
            for i, row in self.annotations.iterrows():
                if row["Feature"]=="CDS":
                    
                    if "ORF_start" in self.annotations.columns and "ORF_end" in self.annotations.columns:
                        kk1,kk2 = "ORF_start","ORF_end"
                    elif "Start" in self.annotations.columns and "Stop" in self.annotations.columns:
                        kk1,kk2 = "Start","Stop"
                        
                    if row["Strand"] in ["F","+"]:
                        gene_dict[range(int(row[kk1]),int(row[kk2]))] = row.Locus
                        feature_dict[range(int(row[kk1]),int(row[kk2]))] = "Coding"
                        gene_dict[range(int(row[kk1]-300),int(row[kk2]))] = "%s-UTR"%row.Locus
                        feature_dict[range(int(row[kk1]-300),int(row[kk2]))] = "Non-coding"
                        
                    elif row["Strand"]=="RC":
                        gene_dict[range(int(row[kk2]),int(row[kk1]))] = row.Locus
                        feature_dict[range(int(row[kk2]),int(row[kk1]))] = "Coding"
                        gene_dict[range(int(row[kk2]+300),int(row[kk1]))] = "%s-UTR"%row.Locus
                        feature_dict[range(int(row[kk2]+300),int(row[kk1]))] = "Non-coding"
                        
                    try:
                        self.get_syn_ns(int(row["ORF_start"]),int(row["ORF_end"]))
                    except:
                        pass
                    
                else:
                    if "FROM" in self.annotations.columns and "TO" in self.annotations.columns:
                        kk1,kk2 = "FROM","TO"
                    elif "Start" in self.annotations.columns and "Stop" in self.annotations.columns:
                        kk1,kk2 = "Start","Stop"
                        
                    gene_dict[range(row[kk1],row[kk2])] = row.Locus
                    feature_dict[range(row[kk1],row[kk2])] = "Non-coding"
                
            self.annotations_dict = RangeDict(gene_dict)
            self.feature_dict = RangeDict(feature_dict)

            with open(self.annotations_pickled_dict,"wb") as file1:
                pickle.dump(self.annotations_dict,file1,protocol=pickle.HIGHEST_PROTOCOL)

            with open(self.features_pickled_dict,"wb") as file2:
                pickle.dump(self.feature_dict,file2,protocol=pickle.HIGHEST_PROTOCOL)

    def annotate_variants(self):

        self.variants_annotated = self.out_directory+"/"+self.sample_name+"_combined_variants_annotated.csv"
        if os.path.isfile(self.variants_annotated):
            pass
        else:
            variants = pd.read_csv(self.variants_combined,index_col=0)
            if self.annotations_dict is None:
                self.create_annotations_dict()
            else:
                pass
            for i2,row2 in variants.iterrows():
                try:
                    annotate = self.annotations_dict[i2]
                    feature = self.feature_dict[i2]
                    variants.at[i2,"Locus"] = annotate
                    variants.at[i2,"Feature"] = feature
                except:
                    pass
                if i2 in self.syn_pos:
                    variants.at[i2,"Codon"] = "S"
                elif i2 in self.ns_pos:
                    variants.at[i2,"Codon"] = "NS"
                else:
                    pass

            variants.to_csv(self.variants_annotated)
    
def identify_ORFs_in_CDS(ref,annotations):
    fasta = SeqIO.read(ref,"fasta")
    ann_pd = pd.read_csv(annotations,delimiter="\t",index_col=0)
    for n,row in ann_pd.iterrows():
        if row.Feature =="CDS":
            F,T = row[["FROM","TO"]].values
            o,l,s,p = ORFpy.extract_orfs(fasta.seq[F:T],F,T)
            try:
                max_l = max(l)
                max_l_i = l.index(max_l)
                ann_pd.at[n,"ORF_start"] = int(p[max_l_i])
                ann_pd.at[n,"ORF_end"] = int(p[max_l_i]+max_l)
                ann_pd.at[n,"Strand"]= s[max_l_i]
            except:
                pass
    return ann_pd

def run_all(sample_folder,pd,ref_fasta,ref_text):
    for reads_folder in os.listdir(sample_folder):
        rf = os.path.join(sample_folder, reads_folder)
        if os.path.isdir(rf):
            reads = [os.path.join(rf,i) for i in os.listdir(rf) if i != '.DS_Store']
            sample_name = reads_folder
            print("Beginning variant calling for %s"%sample_name)
            try:
                a= aligner(ref_fasta,ref_text,sample_name)
                a.parent_directory = pd
                a.reads = reads
                a.check_samples()
                a.create_out_directory()
                a.create_index_bowtie()
                a.create_index_bwa()
                a.align_reads_bowtie()
                print(a.bowtie_coverage)
                a.align_reads_bwa()
                print(a.bwa_coverage)
                a.merge_alignments()
                a.call_variants_bcftools()
                a.call_variants_freebayes()
                a.f_and_c_variants()
                a.annotate_variants()
            except:
                pass
            print("Finished variant calling for %s"%sample_name)
        
        
sf0 = "../../Project_specific_files/RIF_R_evol_project/Evolution_of_tolerance/Sequencing/seq_files"
parent_dir0 = "../../Project_specific_files/RIF_R_evol_project/Evolution_of_tolerance/Sequencing/vcf_files2"
ref0 = "references_and_indexes/MSM_MC2_155.fasta"
ann0 = "references_and_indexes/MSM_MC2_155.txt"

sf1 = "../../Project_specific_files/RIF_R_evol_project/Evolutionary_trajectories/Genomic_DNA_sequencing/MSM_RIF_resistants_GenSeq"
parent_dir1 = "../../Project_specific_files/RIF_R_evol_project/Evolutionary_trajectories/Genomic_DNA_sequencing/MSM_RIF_resistants_VCFs/SNAPPER_calls"

sf2 = "../../Project_specific_files/DRonA_MLSynergy/SA161/Genome_analysis/Sequence_files"
parent_dir2 = "../../Project_specific_files/DRonA_MLSynergy/SA161/Genome_analysis/Variants"
ref2 = "../../Project_specific_files/DRonA_MLSynergy/SA161/Genome_analysis/References/Mtb_H37Rv_genome_v4.fasta"
ann2 = "../../Project_specific_files/DRonA_MLSynergy/SA161/Genome_analysis/References/Mtb_H37Rv_txt_v4.txt"


sf3 = "../../Project_specific_files/RIF_R_evol_project/Clinical_strains_of_Mtb/LLR_Shea/TB_fasta_220311"
parent_dir3 = "../../Project_specific_files/RIF_R_evol_project/Clinical_strains_of_Mtb/LLR_Shea/Variants"
ref3 = "../../Project_specific_files/DRonA_MLSynergy/SA161/Genome_analysis/References/Mtb_H37Rv_genome_v4.fasta"
ann3 = "../../Project_specific_files/DRonA_MLSynergy/SA161/Genome_analysis/References/Mtb_H37Rv_txt_v4.txt"
