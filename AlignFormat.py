import warnings
import numpy as np
import click
import os
import sys
warnings.filterwarnings("ignore")

from Bio import SeqIO, AlignIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


NUC_IND = {'A': 1, 'T': 2, 'C': 3, 'G': 4}


def get_geneslist(alignlist):
    """Function that return the shared list of genes by several genomes"""
    seqlist = set([])
    for al in alignlist:
        seqlist.update([x.id for x in al])
    return seqlist


def concat(alignlist, genename, missing='-', alpha=generic_dna, stopchar='*'):
    """Concatenate alignment into one global alignment. 
    This was mainly copied from coretracker package"""
    seqlist = get_geneslist(alignlist)
    #print(seqlist, len(seqlist))
    genepos = []
    lastpos = 0
    align_dict = {}

    for i, gene in enumerate(genename):  # range(len(alignlist)):
        cur_al = alignlist[i]
        al_len = cur_al.get_alignment_length()
        genepos.append((gene, lastpos, al_len + lastpos))
        lastpos += al_len
        spec_dict = dict((x.id, x) for x in cur_al)
        desc_dict = dict((x.id, x.description) for x in cur_al)
        for spec in seqlist:
            if(spec_dict.get(spec)):
                adding = spec_dict[spec]
            else:
                adding = SeqRecord(Seq(missing * al_len, alpha), id=spec,
                                   name=spec, description=desc_dict.get(spec, ''))
            try:
                align_dict[spec] += adding
            except (KeyError, ValueError):
                align_dict[spec] = adding
    return align_dict, genepos


def nuc_one_shot(x):
    """Utility function that encode each nucleotide as a one-hot vector
    Gap and undefined nucleotides are grouped as a fifth category"""
    val = np.zeros(5)
    val[NUC_IND.get(x, len(val))-1] = 1  # set unknown nuc to last position
    return val


def seq_one_shot(seqrow):
    """Utility function that encode the sequences of an alignment"""
    return np.array([nuc_one_shot(residu) for residu in seqrow])


def save_metadata(gpos, seqpos, outfile):
    """Utility function that save the position of each gene and genome
    inside the numpy matrix"""
    with open(outfile, 'w') as OUT:
        OUT.write('>SeqPos\n')
        OUT.write('#Genome\tRowNumber\n')
        for i, seq in enumerate(seqpos):
            OUT.write('{}\t{}\n'.format(seq, i+1))
        OUT.write('\n>GenePos\n')
        OUT.write('#Gene\tStart\tEnd\n')
        for gene in gpos:
            OUT.write('{}\t{}\t{}\n'.format(*gene))


@click.command(help="Construct and save the representation of the alignment as a matrix")
@click.argument('alignment', nargs=-1, type=click.Path(exists=True))
@click.option('--out', help='Outfile name in which the numpy matrix will be saved')
@click.option('--is_prot', is_flag=True, help='Reference genetic code to use for all genomes.')
def al_to_mat(alignment, out, is_prot=False):
    """Format an alignment and convert it into a one-hot encoding
    that can be used in the algorithms described in the document
    :param alignment: list of alignment files (codon alignment)
    :param out: basename of the output file in which the numpy matrix will be saved
    :return: None
    """
    alphabet = generic_protein if is_prot else generic_dna
    if not isinstance(alignment, (list, tuple)):
        alignment = [alignment]
    alignment = list(alignment)
    genename = []
    for i, al in enumerate(alignment):
        alignment[i] = AlignIO.read(al, format='fasta')
        genename.append(os.path.basename(al).split('.')[0])
    f_align, gpos = concat(alignment, genename, alpha=alphabet)
    f_align = MultipleSeqAlignment(f_align.values())
    seq_position = [seq.name for seq in f_align]
    # reconvert dict to multseq
    seq_rep = np.apply_along_axis(
        seq_one_shot, axis=1, arr=np.asarray(f_align))  # hope axis 1 is correct
    np.save(out, seq_rep, allow_pickle=True)
    save_metadata(gpos, seq_position, os.path.expanduser(out+".info"))


def generate_ed_position(alignment, nuc2sub='C'):
    """Given an alignment generate a matrix containing 
if __name__ == '__main__':
    al_to_mat()
