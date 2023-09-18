#! /usr/bin/env python3
import argparse
import subprocess
from collections import Counter
import pandas as pd
import os

VirBot_path = str(os.path.dirname(os.path.abspath(__file__)))

global positive_cluster
positive_cluster = []

# argument parser
def virbot_cmd():
    parser = argparse.ArgumentParser(description="ARGUMENTS")
    parser.add_argument('--input', type=str, help="The input contig file.")
    parser.add_argument('--output', default="VB_result", type=str, help="The output directory.")
    parser.add_argument('--sen', action='store_true', help="Run the sensitive mode of VirBot.")
    parser.add_argument('--taxa', default="TOP", help="The mode of VirBot's taxanomic module (TOP(default)/LCA)")
    parser.add_argument('--threads', default="8", help="The threads number run for HMMER and DIAMOND")
    args = parser.parse_args()
    return args

def read_thresholding():
    filename = VirBot_path + "/ref/VirBot_hmm_threshold.txt"
    threshold = {}
    with open(filename, 'r') as f:
        for line in f:
            t = line.strip().split()
            threshold[int(t[0])] = float(t[1])
    return threshold
    
def read_hmmtaxa():
    filename = VirBot_path + "/ref/VirBot_hmm_taxa_full.txt"
    hmm_taxa = {}
    with open(filename, 'r') as f:
        for line in f:
            t = line.strip().split('\t')
            hmm_taxa[int(t[0])] = t[1]
    return hmm_taxa

def read_rv_acc():
    filename = VirBot_path + "/ref/VirBot_RNAvirus_acc.txt"
    db_rc_acc = set()
    with open(filename, 'r') as f:
        for line in f:
            db_rc_acc.add(line.strip())
    return db_rc_acc
    
def longest_substring(sa,sb):
    for i in range(min(len(sa),len(sb))):
        if sa[i]!=sb[i]:
         return sa[0:i]
    return sa[0:i]

class contig:
    def __init__(self,fullname):
        self.fullname = fullname
        self.name = fullname.strip('>').split()[0]
        self.seq = None
        self.proteins = {}
        self.rnaviralness = 0
        self.taxa = None

    def add_protein(self,prot_name,prot):
        self.proteins[prot_name] = prot

    def calculate_rnaviralness(self):
#        t, tmp_top_score = 0, 0
        t = 0
        for key, value in  self.proteins.items():
            if value.rnavaralness:
                t += 1
#                if value.potential_taxa and tmp_top_score < value.score:
#                    self.taxa = value.potential_taxa
#                    tmp_top_score = value.score
        if len(self.proteins):
            self.rnaviralness = t / len(self.proteins)
            
    def taxa_assignment(self, taxa_mode = "TOP"):
        if taxa_mode == "TOP":
            tmp_top_score = 0
            for key, value in  self.proteins.items():
                if value.rnavaralness:
                    if value.potential_taxa and tmp_top_score < value.score:
                        self.taxa = value.potential_taxa
                        tmp_top_score = value.score
        elif taxa_mode == "LCA":
            for key, value in self.proteins.items():
                if value.rnavaralness and value.potential_taxa:
                    if self.taxa == None:
                        self.taxa = value.potential_taxa
                    else:
                        self.taxa = longest_substring(self.taxa,value.potential_taxa)

#####################################################################################
class protein:
    db_threshold = read_thresholding()
    db_hmmtaxa = read_hmmtaxa()
    db_rv_acc = read_rv_acc()

    def __init__(self,fullname):
        self.fullname = fullname
        self.contig_name = self._set_contig_name(fullname)
        self.seq=''
        self.potential_taxa = None
        self.best_hit = 0
        self.e_value = 1
        self.score = 0
        self.rnavaralness = 0
        self.diamond = None

    def _set_contig_name(self, fullname):
        t = fullname.strip('>').split()[0]
        t = t.rsplit('_',1)[0]
        return t

    def _set_protein_name(self,fullname):
        t=fullname.split('#')[0][1:-1]
        return t

    def filter2neg(self):
        # only keep the best hit; If it is not pos, then block other hit
        # self.best_hit = 0
        self.e_value = 1
        self.rnavaralness = 0
        self.score = 0

    def parse_match_search(self, result):
        t = result.split()
        hit_clustet_index, e_value, score  \
            = int(t[2].strip("cluster_")), float(t[4]), float(t[5])

        if score > self.score:
            self.score, self.e_value, self.best_hit \
                = score, e_value, hit_clustet_index

    def rnavaralness_for_search(self):
        if self.best_hit and self.score >= protein.db_threshold[self.best_hit] and self.e_value< 1e-3:
            self.rnavaralness = 1
            positive_cluster.append(self.best_hit)
            self.potential_taxa = protein.db_hmmtaxa[self.best_hit]
            print(self.fullname.split('#')[0][1:-1],
                  '\tcluster_%d_(%.1f)' % (self.best_hit, protein.db_threshold[self.best_hit]),
                  '\t', self.score, '\t', self.e_value,
                  '\t', self.potential_taxa)

    def parse_match_diamond_blastx(self, result):
        t = result.split()
        hit_acc = t[1]
        if hit_acc in protein.db_rv_acc:
            self.rnavaralness, self.diamond_acc = \
                1, hit_acc
            print(self.fullname.split('#')[0][1:-1], '\t', hit_acc,
                  '\t', t[-1], '\t', t[-2])

#####################################################################################
def predict(output_dir, temp_dir,
            file_input, file_output,
            sen=False, taxa_mode = "TOP"):
    """
    :param file_input: input contigs
    :param file_output: positive contigs predicted by our strategy
    :return: None
    """
    def parse_protein(filename):
        '''
        :param filename: predicted protein file
        :return: a dict recording all the proteins
        '''
        proteins={}
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('>'):
                    prot_name = line.strip('>').split()[0]
                    proteins[prot_name]=protein(line)
        return proteins

    def parse_hmmsearch(filename,proteins):
        '''
        :param filename: hmmsearch result file
        :param proteins: the tmp protein dict
        :return: the tmp protein dict with hmmsearch result
        '''
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    t = line.strip().split()
                    prot_name = t[0]
                    if prot_name in proteins:
                        proteins[prot_name].parse_match_search(line)

        print("Protein_acc\tHMM_name_(corresponding_threshold)\tBit_score\tE-value\tTaxa")
        for prot_name, protein in proteins.items():
            protein.rnavaralness_for_search()

        return proteins

    def parse_diamond_blastx(filename, proteins):
        '''
        :param filename: the DIAMOND result diamond
        :param proteins: the tmp protein dict with prediction result
        :return: the tmp protein dict with sensitive prediction result
        '''
        print("Protein_acc\tReference_acc\tBit_score\tE-value")
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    t = line.strip().split()
                    prot_name = t[0]
                    if prot_name in proteins:
                        proteins[prot_name].parse_match_diamond_blastx(line)
        return proteins

    def parse_contig(filename,proteins):
        '''
        :param filename: the input file
        :param proteins: the tmp protein dict
        :return: the tmp contigs dict that containing description and encoded proteins
        '''
        contigs={}
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_name = line.strip('>').split()[0]
                    contigs[contig_name] = contig(line)

        for prot_name, prot in proteins.items():
            if prot.contig_name in contigs:
                contigs[prot.contig_name].add_protein(prot_name,prot)

        for c_name, c in contigs.items():
            c.calculate_rnaviralness()

        return contigs

    def output_rnaviralness(filename, contigs, taxa_mode):
        """
        :param filename: the input file
        :param contigs: the tmp contigs dict that containing description and encoded proteins
        :return: the tmp contigs dict that containing the accession and protein prediction result
        """
        i = 0
        print('Contig_acc\tRNA-viral_gene_content\tEncoded_proteins_num\tLikely_taxa')
        positive_contigs = {}
        # Sum the proteins into the contig, and recode the acc.
        for c_name,c in contigs.items():
            c.taxa_assignment(taxa_mode)
            # viral gene content cutoff: 1/16
            if c.rnaviralness >= 0.0625:
                if c.taxa == None:
                    c.taxa = "Riboviria"
                print(c_name,'\t',round(c.rnaviralness,2),'\t',len(c.proteins),'\t',c.taxa)
                positive_contigs[c_name] = c
            if c.proteins.__len__()>0:
                i+=1

        print("total num of contigs:",contigs.__len__())
        print("num of contigs containing gene(s):",i)
        print("num of contigs lacking gene(s)",contigs.__len__()-i)
        print("num of positive_contigs:",positive_contigs.__len__())
        # print("Counter of cluster hit",Counter(positive_cluster))

        # Retrieve the DNA sequences of positive contigs.
        with open(filename,'r') as f:

            significant =False
            seq = ''
            contig_name = None

            for line in f:
                if line.startswith('>'):
                    if significant:
                        positive_contigs[contig_name].seq= seq

                    significant = False
                    contig_name = line.strip('>').split()[0]
                    if contig_name in positive_contigs:
                        significant = True
                        seq=''
                elif significant:
                    seq += line
            #do not ignore the last contig
            if significant:
                positive_contigs[contig_name].seq = seq
        return positive_contigs

    def write_positive_file(output_dir, file_output, positive_contigs):
        """
        :param filename: the output filename
        :param positive_contigs: all the predicted positive contigs
        :return: None
        """
        positive_contig_output = []
        with open(output_dir + file_output,'w') as f:
            for c_name,c in positive_contigs.items():
                f.write(c.fullname)
                f.write(c.seq)
                positive_contig_output.append([c_name, round(c.rnaviralness,2), len(c.proteins), c.taxa])

        df = pd.DataFrame(positive_contig_output, columns=['Contig_acc','RNA-viral_gene_content','Encoded_proteins_num','Likely_taxa'])
        df.to_csv(output_dir + "pos_contig_score.csv",index=False)

#####################################################################################
    # store all the predicted proteins into a ditc {prot_name: proteins_information} without the seq.
    proteins = parse_protein(temp_dir + "protein.faa")
    print("Num of proteins",proteins.__len__())

    # parse the hmmsearch result for the proteins, ditc {prot_name: proteins_information}.
    proteins = parse_hmmsearch(temp_dir + "VB_hmmer.out", proteins)
    print("Parsing of protein HMM-match result finished")

    # parse the diamond blasp result for the proteins, ditc {prot_name: proteins_information}.
    if sen:
       proteins = parse_diamond_blastx(temp_dir + "VB_diamond.out", proteins)
       print("Parsing of protein DIAMOND-match result finished")

    # record all the contigs from input return, ditc {contig_name: contigs_information}.
    contigs = parse_contig(file_input, proteins)

    # summarize the protein result into contigs, and retrieve the identified sequences.
    # ditc {positive_contig_name: contigs_information}.
    positive_contigs = output_rnaviralness(file_input, contigs, taxa_mode)

    # write the output file
    write_positive_file(output_dir, file_output, positive_contigs)

#####################################################################################
if __name__ == "__main__":

    args = virbot_cmd()
    
    if args.taxa != "TOP" and args.taxa != "LCA":
        raise Exception("The taxonmic mode should be \"TOP\" or \"LCA\".")

    output_dir = args.output
    if os.path.exists(output_dir):
        raise Exception("The output directory already exists.")
    else:
        os.makedirs(output_dir)
    temp_dir = f"{output_dir}/tmp"
    os.makedirs(temp_dir)

    FNULL = open(os.devnull, 'w')

    print(f"Input contig: {args.input}")
    print(f"Output directory: {output_dir}")
    if args.sen:
        print(f"Mode: sensitive mode")
    else:
        print(f"Mode: default mode")

    # run Prodigal
    print("Predicting the encoded proteins...")
    subprocess.run(f"prodigal -i {args.input} -a {temp_dir}/protein.faa -p meta", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    print("Proteins prediction finished.")

    # run HMMER
    print("Scanning the protein by hmmsearch...")
    subprocess.run(f"hmmsearch --tblout {temp_dir}/VB_hmmer.out --noali -E 0.001 --cpu {args.threads} {VirBot_path}/ref/VirBot.hmm {temp_dir}/protein.faa", shell=True, stdout=FNULL)
    print("HMMER finshed.")

    # run DIAMOND (in sensitive mode)
    if args.sen:
        print("Scanning the protein by DIAMOND...")
        subprocess.run(f"diamond blastp --db {VirBot_path}/ref/VirBot.dmnd --query {temp_dir}/protein.faa "
                       f"--outfmt 6 --max-target-seqs 1 --threads {args.threads} "
                       f"--evalue 1e-5 --out {temp_dir}/VB_diamond.out", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        print("DIAMOND finshed.")

    # predict using VirBot
    predict(output_dir=f"{output_dir}/",
            temp_dir=f"{temp_dir}/",
            file_input = args.input,
            file_output = "output.vb.fasta",
            sen=args.sen, taxa_mode = args.taxa)

