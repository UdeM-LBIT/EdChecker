import click
import os
import pyvolve
import sys
import numpy as np
from ete3 import Tree
from collections import defaultdict as ddict
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord

GTABLE = CodonTable.unambiguous_dna_by_id[1]


def isProb(*args):
    try:
        val = float(args[0])
        if not 0 <= val <= 1.0:
            raise ValueError()
        return val
    except:
        raise ValueError("{} should be a float between 0-1".format(args[0]))


def mkdir(dirname):
    """Clean a directory"""
    os.makedirs(dirname, exist_ok=True)
    return dirname


def add_height_to_tree(tree):
    for node in tree.traverse():
        node.add_features(height=node.get_distance(tree, topology_only=True))
    return tree


class Mutex(click.Option):
    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop("not_required_if", [])
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help", "") + "\nOption is mutually exclusive with " +
                          ", ".join(self.not_required_if) + ".").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '" + str(self.name) + "' is mutually exclusive with " + str(mutex_opt) + ".")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


@click.command(help="Simulation of a dataset with C2U RNA editing")
@click.option('--gnumber', type=click.INT, cls=Mutex,  not_required_if=["glist", "gtree"], help='Number of genomes in the dataset. This is provided for convenience, but you should use --glist or --gtree as they give more freedom.')
@click.option('--glist', type=click.Path(exists=True), cls=Mutex,  not_required_if=["gnumber", "gtree"], help='List of leaf (genome) names. Start a genome with \'-\' or \'_\'  if you don\'t want RNA editing in it. Comment lines start with \'#\'')
@click.option('--gtree', type=click.Path(exists=True), cls=Mutex,  not_required_if=["gnumber", "glist"], help='A phylogenetic tree that depicts the evolutionary relationship between genomes. Branch lengths should be present ! Start a genome with \'_\' if you do not want RNA editing in it.')
@click.option('--gsize', type=click.INT, help='Number of genes per genome (see it as the number of simulations)')
@click.option('--edprob', default=0.1, type=isProb, help='RNA editing probability at a random site, per genome. The number of edited site is a taken from a normal distribution ')
@click.option('--glen_range', default=(200, 500), nargs=2, type=click.Tuple([int, int]), help='Min and Max length (as number of AA) of each coding gene of genomes')
@click.option('--dnds', default=(0.5, 1), nargs=2, type=click.Tuple([float, float]), help='Rate of synonymous and nonsynonymous changes (dN, dS)')
@click.option('--tau', type=click.FLOAT, help='Transition / Transversion rate')
@click.option('--delrate', type=isProb, help='Sequence random deletion rate for codons (Pyvolve does not allow indel, so this randomly choose some position of the alignment then delete the amino acid.')
@click.option('--from_al', type=click.Path(exists=True), help='Compute codon frequencies from existing alignment')
@click.option('--protlike', is_flag=True, help='Automatically add the ATG methionine codon at the start of each sequences.')
@click.option('--no_syn', is_flag=True, help='Ignore RNA editing positions that results just in synonymous codons. This will make optimization-like methods simpler')
@click.option('--sub_rate', type=click.FLOAT, default=1, help='Substitution rates. Branch length will be rescale by this factor. Increase or decrease to change the overall mutation rate.')
@click.option('--min_cons', type=isProb, default=0.0, help='Minimum column conservation (0-1 as in percent of maximum conservation) for a column to be considered for RNA editing.')

@click.option('-wd', '--outdir',  default="./dataset/simul", help='Output directory where to save the dataset and all results')
def cli(gnumber, glist, gtree, edprob, gsize, glen_range, dnds, tau=None, delrate=0.0, from_al=None, protlike=False, no_syn=False, sub_rate=1.0, min_cons=0.0, outdir=""):
    """Extract genome content based on a list of species """
    gleaf = []
    no_edit = []
    tree = None
    if gnumber:
        gleaf = ['Genome_{}'.format(i) for i in range(1, gnumber+1)]
    elif glist:
        with open(glist) as G:
            for line in Glist:
                line = line.strip()
                if line and not line.startswith('#'):
                    gleaf.append(line.strip('-_'))
                    if line.startswith('-') or line.startswith('_'):
                        no_edit.append(line.strip('-_'))
    elif gtree:
        tree = Tree(gtree)
        gleaf = tree.get_leaf_names()
        no_edit = [x.strip('_') for x in gleaf if x.startswith('_')]
        for node in tree:
            node.name = node.name.strip('_')

    else:
        raise NotImplementedError(
            "One of --gnumber, --glist and --gtree is needed !")

    if not tree:
        tree = Tree()
        tree.populate(len(gleaf), names_library=gleaf, random_branches=True)

    param_list = {"alpha": dnds[1], "beta": dnds[0]}
    if tau:
        param_list.update({"kappa": tau})

    if from_al:  # read codons frequencies from an existing alignment
        f = pyvolve.ReadFrequencies("codon", file=from_al)
        param_list.update({'state_freqs': f.compute_frequencies()})

    #print(tree.get_ascii(show_internal=True, attributes=['name', 'dist']))
    phylogeny = pyvolve.read_tree(tree=tree.write(format=5), scale_tree=sub_rate)
    codon_model = pyvolve.Model("codon", param_list)#, neutral_scaling=True)
    sequences = []
    edited_sequences = []
    truth_tables = []
    # add height to tree
    tree = add_height_to_tree(tree)

    for i in range(gsize):
        # gene length is given from an uniform distribution
        alen = np.random.randint(glen_range[0], glen_range[1])*3
        seq = simulate_genomes(codon_model, phylogeny, alen, outdir, i+1)
        if delrate:
            seq = random_deletion(seq, tree, alen//3, delrate)
        if protlike:
            for k in seq:
                seq[k] = 'ATG'+seq[k]
        sequences.append(seq)
        edited_seq, truth_table = CtoUsimulate(
            seq, tree, no_edit, edprob, no_syn=no_syn, min_cons=min_cons)
        edited_sequences.append(edited_seq)
        truth_tables.append(truth_table)
        save_data(tree, seq, edited_seq, truth_table, outdir, i+1)


def save_data(tree, seq, edseq, ttable, outdir, number):
    path = mkdir(os.path.join(outdir, str(number)))
    with open(os.path.join(path, "rna_{}.fasta".format(number)), 'w') as OUT:
        for k, v in seq.items():
            OUT.write('>{}\n{}\n'.format(k, v))
    with open(os.path.join(path, 'dna_{}.fasta'.format(number)), 'w') as EOUT:
        for k, v in edseq.items():
            EOUT.write('>{}\n{}\n'.format(k, v))
    with open(os.path.join(path, 'dna_{}.truth'.format(number)), 'w') as TVAL:
        for spec, pos in ttable.items():
            if pos:
                TVAL.write('{}: {}\n'.format(
                    spec, ",".join((str(x) for x in pos))))
    tree.write(outfile=os.path.join(outdir, 'tree.nw'))


def scan_for_nuc(seq, nuc='T'):
    seqid, seqval = zip(*seq.items())
    np_align = np.asarray([[snuc for snuc in al] for al in seqval])
    return np_align, np_align == nuc, dict((val, i) for i, val in enumerate(seqid))


def codon2aa(codon):
    if codon == '---':
        return '-'
    elif codon in GTABLE.stop_codons:
        return '*'
    return GTABLE.forward_table.get(codon, 'X')


def translate_seq(seqs):
    stop_codons = GTABLE.stop_codons
    aa_seqs = []
    for seqid, seq in seqs.items():
        prot_seq = [codon2aa(seq[i:i+3]) for i in range(0, len(seq), 3)]
        aa_seqs.append(SeqRecord(Seq("".join(prot_seq), alphabet=generic_protein),id=seqid))
    msa = MultipleSeqAlignment(aa_seqs)
    return msa


def compute_ic_content(alignment):
    align_info = AlignInfo.SummaryInfo(alignment)
    align_info.information_content()
    ic_vector = align_info.ic_vector
    # biopython wrong version hack
    if isinstance(ic_vector, dict):
        ic_vector = np.zeros(len(align_info.ic_vector))
        for (ic_i, ic_v) in align_info.ic_vector.items():
            ic_vector[ic_i] = ic_v
    return list(ic_vector)


def CtoUsimulate(seq, tree, no_edit, edprob, no_syn=False, min_cons=0.0):
    # only one RNA editing per codons
    prot_msa = translate_seq(seq)
    ic_content = compute_ic_content(prot_msa)
    max_ic_content = max(ic_content)
    nodelist, nodeheight = get_nodelist(tree, outside=no_edit)
    seq_align, T_positions, index_getter = scan_for_nuc(seq)
    seq_len = seq_align.shape[1]
    truth_table = ddict(list)
    for pos in range(0, seq_len, 3):
        selected_node = np.random.choice(nodelist, p=nodeheight)
        spec_list = [x.name for x in selected_node if x.name not in no_edit]
        current_col_cons = ic_content[pos//3]
        if max_ic_content*min_cons<=current_col_cons and spec_list and np.random.rand() < edprob:
            #print(current_col_cons, pos)
            spec_array = np.asarray([index_getter[x] for x in spec_list])
            best_codon_pos = np.sum(T_positions[spec_array, pos:pos+3], axis=0)
            best_codon_pos = np.argmax(best_codon_pos)
            for spec in spec_list:
                if seq_align[index_getter[spec], pos+best_codon_pos] == 'T':
                    seq_align[index_getter[spec], pos+best_codon_pos] = 'C'
                    truth_table[spec].append(pos+best_codon_pos)

    new_align = {}
    for spec, ind in index_getter.items():
        if no_syn:
            unwanted = []
            for ed_pos in truth_table[spec]:
                nuc_pos = (ed_pos//3)*3
                codon_ed = seq[spec][nuc_pos: nuc_pos+3]
                codon = "".join(seq_align[ind, nuc_pos:nuc_pos+3])
                aa = GTABLE.forward_table.get(codon_ed)
                if aa and aa == GTABLE.forward_table.get(codon):
                    # PLEASE COMMENT THE NEXT LINE IF YOU DO NOT CARE ABOUT
                    # RESETING THE CODON ITSELF, BUT JUST WANT TO
                    # REMOVE THE POSITION
                    seq_align[ind, ed_pos] = 'T'  # reset to original val
                    # reset to prev value
                    unwanted.append(ed_pos)
            truth_table[spec] = [
                x for x in truth_table[spec] if x not in unwanted]
        new_align[spec] = "".join(seq_align[ind, :])
        # we do not want to keep éééé

    return new_align, truth_table


def get_nodelist(tree, outside=[]):
    # if outside is provided
    # just restricted to tree root to LCA of accepted genomes
    best_lca = tree
    if outside:
        best_lca = tree.get_common_ancestor([x for x in tree if x.name not in outside])
    all_height = sum([x.height for x in best_lca.get_descendants()])
    nodelist, nodeheight = zip(
        *[(x, x.height*1.0/all_height) for x in best_lca.get_descendants()])
    return nodelist, nodeheight


def random_deletion(seq, tree, asize, delrate=0.05):
    nodelist, nodeheight = get_nodelist(tree)
    for k, v in seq.items():
        seq[k] = [x for x in v]
    for pos in range(0, asize, 3):
        prob = np.random.rand()
        if prob <= delrate:
            # we should delete column
            seq_with_del = np.random.choice(nodelist, p=nodeheight)
            for l in seq_with_del.get_leaf_names():
                seq[l][pos:pos+3] = '---'
    for k, v in seq.items():
        seq[k] = "".join(v)
    return seq


def simulate_genomes(model, tree, asize, outdir, number):
    path = mkdir(os.path.join(outdir, str(number)))
    partition = pyvolve.Partition(models=model, size=asize)
    evolver = pyvolve.Evolver(tree=tree, partitions=partition)
    evolver(seqfile=None,  # ,
            ratefile=os.path.join(path, "rate_{}.fasta".format(number)),
            infofile=None)
    return evolver.get_sequences()


if __name__ == '__main__':
    cli()
