import click
from ete3 import Tree
import pyvolve
import os
import numpy as np
import sys

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
        self.not_required_if:list = kwargs.pop("not_required_if", [])
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help", "") + "\nOption is mutually exclusive with " + ", ".join(self.not_required_if) + ".").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt:bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError("Illegal usage: '" + str(self.name) + "' is mutually exclusive with " + str(mutex_opt) + ".")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)

@click.command(help="Simulation of a dataset with C2U RNA editing")
@click.option('--gnumber', type=click.INT, cls=Mutex,  not_required_if=["glist", "gtree"], help='Number of genome to use')
@click.option('--glist', type=click.Path(exists=True), cls=Mutex,  not_required_if=["gnumber", "gtree"], help='List of leaf (genome) names. Start a genome with \'-\' if you don\'t want RNA editing in it.' )
@click.option('--gtree', type=click.Path(exists=True), cls=Mutex,  not_required_if=["gnumber", "glist"], help='A phylogenetic tree that depict the evolutionary relationship between genomes')
@click.option('--gsize', type=click.INT, help='Number of genes per genome')

@click.option('--edprob', default=0.1, type=isProb, help='RNA editing probability at a random site, per genome. The number of edited site is a taken from a normal distribution ')
@click.option('--glen_range', default=(200, 500), nargs=2, type=click.Tuple([int, int]), help='Min and Max length (as number of AA) of each coding gene of genomes')
@click.option('--dnds', default=(0.5, 1), nargs=2, type=click.Tuple([float, float]), help='Rate of synonymous and nonsynonymous changes (dN, dS)')
@click.option('--tau', type=click.FLOAT, help='Transition / Transversion rate')
@click.option('--delrate', type=isProb, help='Sequence random deletion rate for codons (Pyvolve does not allow indel, so this randomly choose some position of the alignment then delete the amino acid.')
@click.option('--from_al', type=click.Path(exists=True), help='Compute codon frequencies from existing alignment')
@click.option('-wd', '--outdir',  default="./dataset/simul", help='Output directory where to save the dataset and all results')
def cli(gnumber, glist, gtree, edprob, gsize, glen_range, dnds, tau=None, delrate=0.0, from_al=None, outdir=""):
    """Extract genome content based on a list of species """
    print(delrate)
    gleaf = []
    no_edit = []
    tree = None
    if gnumber:
        gleaf = ['Genome_{}'.format(i) for i in range(1, gnumber+1)]
    elif glist:
        with open(glist) as G:
            for line in Glist:
                line = line.strip()
                if not line.startswith('#'):
                    gleaf.append(line)
                    if line.startswith('-') or line.startswith('_'):
                        no_edit.append(line.strip('-_'))
    elif gtree:
        tree = Tree(gtree)
        gleaf = gtree.get_leaf_names()
        no_edit = [x for x in gleaf if x.startswith('_')]

    else:
        raise NotImplementedError("One of --gnumber, --glist and --gtree is needed !")

    if not tree:
        tree = Tree()
        tree.populate(len(gleaf), names_library=gleaf, random_branches=True)

    param_list = {"alpha":dnds[1], "beta":dnds[0]}
    if tau:
        param_list.update({"kappa":tau})

    if from_al: # read codons frequencies from an existing alignment
        f = pyvolve.ReadFrequencies("codons", file=from_al)
        param_list.update({'state_freqs':f.compute_frequencies()})

    #print(tree.get_ascii(show_internal=True, attributes=['name', 'dist']))
    print(tree.write(format=5))
    phylogeny = pyvolve.read_tree(tree=tree.write(format=5))
    codon_model = pyvolve.Model("codon", param_list, neutral_scaling=True)
    sequences = []
    edited_sequences = []
    # add height to tree
    tree = add_height_to_tree(tree)

    for i in range(gsize):
        alen = np.random.randint(glen_range[0], glen_range[1])*3 # gene length is given from an uniform distribution
        seq = simulate_genomes(codon_model, phylogeny, alen, outdir, i)
        if delrate:
            seq = random_deletion(seq, tree, alen//3, delrate)
        sequences.append(seq)
        edited_seq = CtoUsimulate(seq, tree, no_edit, edprob)
        edited_sequences.append(edited_seq)


def CtoUsimulate(seq, tree, no_edit, edprob):
    print(seq)
    sys.exit(0)
    C_position = scan_for_Cytosine(seq)


def random_deletion(seq, tree, asize, delrate=0.05):
    all_height = sum([x.height for x in tree.traverse()])
    nodelist, nodeheight = zip(*[(x, x.height*1.0/all_height) for x in tree.traverse()])
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
    evolver(seqfile=os.path.join(path, "seq{}.fasta".format(number)), 
            ratefile=os.path.join(path, "rate{}.fasta".format(number)),
            infofile=None)
    return evolver.get_sequences()


if __name__ == '__main__':
    cli()