#/usr/bin/env python
import warnings
import click
import logging
import sys
import os
import subprocess
import re
import sqlite3
warnings.filterwarnings("ignore")  # dangerous

from collections import defaultdict as ddict
from Bio import SeqIO, AlignIO, codonalign
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq
from Bio.Alphabet import IUPAC, generic_nucleotide, generic_protein
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.codonalign.codonalphabet import (
    default_codon_alphabet, get_codon_alphabet)
from Bio.Data import CodonTable

logging.basicConfig(format='DataParser: %(message)s', level=logging.INFO)

EDITING_OK = {'A': ['A'], 'G': ['G'], 'C': [
    'C', 'T'], 'T': ['T'], 'Y': ['T', 'C']}

CtoU = re.compile('C\s*to\s*U', re.IGNORECASE)

MSG = set([])


def mkdir(dirname):
    """Clean a directory"""
    os.makedirs(dirname, exist_ok=True)
    return dirname


def execute_alignment(cmdline, inp, out):
    """Construct alignment command line and execute it """
    prog = 'mafft'
    if 'muscle' in cmdline:
        prog = 'muscle'
        cmdline += " -in %s -out %s" % (inp, out)
    elif 'mafft' in cmdline:
        cmdline += " %s > %s" % (inp, out)
    else:
        raise ValueError(
            "Cannot execute %s. Programme not expected. You can provide your own alignment instead.")
    p = subprocess.Popen(
        cmdline, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    if err:
        logging.debug(err)
    if out:
        pass
    return err


def align(gname, seqlist, msaprog, outdir, alpha=generic_protein, clean_dt=True):
    """Align sequences using muscle of mafft"""
    outdir = mkdir(os.path.expanduser(outdir))
    tmpseq = os.path.join(outdir, "{}.fasta".format(gname))
    align_seq = os.path.join(outdir, "{}.aln".format(gname))
    SeqIO.write(seqlist, open(tmpseq, 'w'), 'fasta')
    execute_alignment(msaprog, tmpseq, align_seq)
    msa = AlignIO.read(align_seq, 'fasta', alphabet=alpha)
    if clean_dt:
        os.remove(tmpseq)
    return msa


def translate_seq(seq, table=1):
    """Useless function to translate sequences
    Only used to test if biopython consider edited site in its translation: Spoiler: NO"""
    ctable = CodonTable.unambiguous_dna_by_id[table]
    stop_codons = ctable.stop_codons
    codon_rep = [seq[i:i+3] for i in range(0, len(seq), 3)]
    prot_seq = []
    for i, cod in enumerate(codon_rep):
        if cod in stop_codons and len(codon_rep) - 1 <= i:
            prot_seq.append('*')
        elif cod == '---':
            prot_seq.append('-')
        else:
            prot_seq.append(ctable.forward_table.get(cod, 'X'))
    return Seq("".join(prot_seq), alphabet=generic_protein)


def load_genbank_file(gbfiles):
    """Get a genbank file"""
    gfile2rec = {}
    for gfile in gbfiles:
        with open(gfile, 'rU') as INFILE:
            rec = SeqIO.read(INFILE, "genbank")
            gfile2rec[gfile] = rec

    return gfile2rec


def load_rev_genes(revfile):
    """Load an alias file that match each gene name to its known synonyms"""
    revgenes = {}
    with open(revfile) as RV:
        for line in RV:
            line = line.strip().split()
            if len(line) > 1:
                for entry in line[1:]:
                    revgenes[entry] = line[0]
    return revgenes


def genematch(gene, glist):
    """Check if gene is found in glist"""
    return any([x.match(gene) for x in glist])


def parse_genelist(glistfile):
    """Parse a list file containing a list of genes (regexp allowed) and return a regexp for each gene"""
    genelist = []
    if glistfile:
        with open(glistfile) as Glist:
            for line in Glist:
                if not line.startswith('#'):
                    genelist.append(re.compile(
                        line.strip().lower().replace('.*', '*').replace('*', '.*')))
    return genelist


def editing_yielder(cod):
    """Iter through all possible edited version of a codon
    We only consider the C --> U editing type"""
    for ednuc0 in EDITING_OK.get(cod[0], [cod[0]]):
        for ednuc1 in EDITING_OK.get(cod[1], [cod[1]]):
            for ednuc2 in EDITING_OK.get(cod[2], [cod[2]]):
                yield ednuc0+ednuc1+ednuc2


def codon_aligner(gname, protalign, genelist, translist, knowneds={}, codemap={}, ignore_syn=False, force_extend=True):
    """DNA alignment codon by codon using protein alignment as template
    Frameshifting is not considered"""
    #aligndict = protalign.to_dict()
    # for seqname, seqrec in aligndict:
    codon_align = ddict(str)
    editing_pos = ddict(list)
    for aseq in protalign:
        glen2alen = {}
        gcode = codemap.get(aseq.name, 1)
        ctable = CodonTable.unambiguous_dna_by_id[gcode]
        codseq = ""

        ungap_len = len(aseq.seq.ungap('-'))
        nuc_seq = genelist[aseq.name]
        if len(nuc_seq) == (ungap_len+1)*3:
            nuc_seq = nuc_seq[:-3]
        trans = translist[aseq.name]
        spec_ed = knowneds.get(aseq.name, [])
        if len(trans) == (ungap_len+1):
            trans = trans[:-1]

        trans_pos = 0
        for apos, val in enumerate(aseq):
            if val == '-':
                codseq += '---'
            else:
                cod = nuc_seq[3*trans_pos:3*trans_pos+3].seq._data
                if not trans[trans_pos] == val and force_extend:  # nothing to say here ?
                    editerator = editing_yielder(cod)
                    while True:
                        try:
                            codmut = next(editerator)
                            if codmut and ctable.forward_table.get(codmut, 'X') == val:
                                editing_pos[aseq.name].extend(
                                    [(npos+trans_pos*3, len(codseq)+npos) for npos, mnuc in enumerate(codmut) if cod[npos] != mnuc])
                                break
                        except StopIteration:
                            break

                glen2alen[trans_pos] = len(codseq)
                codseq += cod
                trans_pos += 1
        # here try to add missing editing position catched by genbank scan or db
        assert len(codseq) == len(aseq)*3, 'Something is wrong'
        for ced in spec_ed:
            ced_cod_pos = ced % 3
            editing_pos[aseq.name].append((ced, glen2alen[ced//3]+ced_cod_pos))

        # new filter on editing_pos: to clean position that do not match with C. Using this and permissive is good
        not_C_matcher = [x for x in editing_pos[aseq.name]
                         if codseq[x[-1]] != 'C']
        editing_pos[aseq.name] = [x for x in sorted(
            set(editing_pos[aseq.name])) if x not in not_C_matcher]

        if ignore_syn:
            rempos = []
            for i, edpos in enumerate(editing_pos[aseq.name]):
                _, edpos = edpos
                cpos = edpos // 3
                old_codon = str(codseq[cpos*3: cpos*3+3])
                new_codon = [x for x in old_codon]
                new_codon[edpos % 3] = 'T'
                new_codon = "".join(new_codon)
                old_aa = ctable.forward_table.get(old_codon, 'X')
                if old_aa == ctable.forward_table.get(new_codon, 'X') and old_aa != 'X':
                    rempos.append(i)
            editing_pos[aseq.name] = [val for i, val in enumerate(
                editing_pos[aseq.name]) if i not in rempos]

        codon_align[aseq.name] = SeqRecord(Seq(
            codseq, alphabet=default_codon_alphabet), id=nuc_seq.id, name=nuc_seq.name, description=nuc_seq.description)
        # if not str(aseq.seq.ungap('-'))==str(trans.seq):
        #     print(editing_pos[aseq.name])
        #     print(aseq.name)
        #     print(ungap_len, len(trans), len(nuc_seq)//3)
        #     assert len(editing_pos[aseq.name]) >0, "At least one position should be edited"

    # codonalign.CodonAlignment(codon_align.values(), alphabet=default_codon_alphabet)
    return editing_pos, codon_align


def get_correct_seq(codon_align, edited_pos, gcode={}):
    """Compute correct sequence from editing position"""
    new_seqs = {}
    for seqid, seqrec in codon_align.items():
        edpos = edited_pos[seqid]
        seqx = [x for x in seqrec.seq._data]
        for pos in edpos:
            seqx[pos[-1]] = 'T'

        new_seqs[seqid] = SeqRecord(translate_seq("".join(seqx), gcode.get(
            seqid, 1)), id=seqrec.id, name=seqrec.name, description=seqrec.description)
    return new_seqs.values()


def at_least_one_ed(edited_pos, pos):
    for s_id, ed_pos_list in edited_pos.items():
        for _, ed_pos in ed_pos_list:
            if ed_pos//3 == pos:
                return ed_pos % 3
    return None


def complete_from_consensus(true_seq_aln, cod_align, edited_pos, gcode={}, only_ed_G=False):
    true_seq_aln = MultipleSeqAlignment(list(true_seq_aln))
    summary_align = AlignInfo.SummaryInfo(true_seq_aln)
    consensus = summary_align.dumb_consensus(threshold=0.5)
    for pos in range(true_seq_aln.get_alignment_length()):
        cons_aa = consensus[pos]
        known_ed_col = at_least_one_ed(edited_pos, pos)
        if known_ed_col is not None:
            for seqrec in true_seq_aln:
                ctable = CodonTable.unambiguous_dna_by_id[gcode.get(
                    seqrec.name, 1)]
                seq_aa = seqrec[pos]
                ed_allowed = (only_ed_G and len(edited_pos.get(
                    seqrec.name, [])) > 0) or not only_ed_G
                nuc_pos = len(str(seqrec[:pos+1].seq).replace("-", ""))
                in_pos = [x[-1]
                          for x in edited_pos[seqrec.name] if x[-1]//3 == pos]
                if in_pos:
                    wg_pos = in_pos[0]
                    # Here we attempt to slightly correct a wrong position
                    wg_cod = [x for x in str(
                        cod_align[seqrec.name][pos*3: pos*3+3].seq)]
                    if not wg_cod[wg_pos % 3] == 'C':
                        edited_pos[seqrec.name] = [
                            x for x in edited_pos[seqrec.name] if x != wg_pos]
                    else:
                        wg_cod[wg_pos % 3] = 'T'
                        wg_aa = ctable.forward_table.get("".join(wg_cod), 'X')
                        if wg_aa != 'X' and wg_aa != cons_aa:
                            edited_pos[seqrec.name] = [
                                x for x in edited_pos[seqrec.name] if x != wg_pos]

                    #print(seqrec.description, pos, pos*3, ed_allowed, known_ed_col, cons_aa, seq_aa, str(cod_align[seqrec.name][pos*3: pos*3+3].seq))

                if ed_allowed and not in_pos:
                    codon = [x for x in str(
                        cod_align[seqrec.name][pos*3: pos*3+3].seq)]
                    cmut = [x for x in codon]
                    cmut[known_ed_col] = 'T'
                    if codon[known_ed_col] == 'C' and ctable.forward_table.get("".join(cmut), 'X') == cons_aa and cons_aa != 'X':
                        edited_pos[seqrec.name].append(
                            (nuc_pos*3+known_ed_col, pos*3 + known_ed_col))
                    else:
                        editerator = editing_yielder(codon)
                        while True:
                            try:
                                codmut = next(editerator)
                                if codmut and ctable.forward_table.get(codmut, 'X') == cons_aa:
                                    edited_pos[seqrec.name].extend(
                                        [(npos+nuc_pos*3, pos*3+npos) for npos, mnuc in enumerate(codmut) if codon[npos] != mnuc])
                                    break
                            except StopIteration:
                                break

    for k, v in edited_pos.items():
        edited_pos[k] = list(sorted(set(v)))
    return edited_pos


def save_data(gname, align, editpos, outdir, trueseq=None):
    """Save data"""
    outdir = mkdir(os.path.expanduser(outdir))
    if trueseq:
        with open(os.path.join(outdir, gname+".aa.aln"), 'w') as AA:
            SeqIO.write(trueseq, AA, 'fasta')
    with open(os.path.join(outdir, gname+".cod.aln"), 'w') as COD:
        SeqIO.write(align, COD, 'fasta')
    with open(os.path.join(outdir, gname+'.truth'), 'w') as TVAL:
        for spec, pos in editpos.items():
            if pos:
                TVAL.write('{}: {}\n'.format(
                    spec, ",".join((str(x[-1]) for x in pos))))


def triple_format(gene2spec, prot2spec, trans2spec, glist=[]):
    """Utility function to quickly print data related to a gene"""
    for k, v in gene2spec.items():
        if not glist or k in glist:
            print('>>'+k)
            for kk, vv in v.items():
                print('++'+kk)
                print('1-Gene'+"\n"+vv.format('fasta'))
                print('1-Prot'+"\n"+prot2spec[k][kk].format('fasta'))
                print('1-Trans'+"\n"+trans2spec[k][kk].format('fasta'))


def collect_edpos_from_seqrec(seqrec):
    """Collect the list of annoted editing positions from a seqrec file. Support only rna editing by nucleotide conversion
    :param seqrec: A seqrecord file
    :return: A dictionnary containing all edited positions for each of the genes in the seqrec
    :rtype: dict
    """
    misc_features = [x for x in seqrec.features if x.type ==
                     "misc_feature" and "editing" in "".join(x.qualifiers.get("note", [])).lower()]
    R_edit_features = [x for x in seqrec.features if x.type ==
                       "RNA_editing" and "substitution" in "".join(x.qualifiers.get("type", [])).lower()]
    edpos = []
    for mf in misc_features:
        affected_gene = mf.qualifiers.get("gene", [None])[0]
        edtype = mf.qualifiers.get("note", [""])[0].replace(
            "RNA editing", "").strip()
        if CtoU.match(edtype):
            edpos.append((affected_gene, mf.location, False))
    if R_edit_features or True:
        for mf in R_edit_features:
            affected_gene = mf.qualifiers.get("gene", [None])[0]
            if mf.qualifiers.get('replace', [""])[0].strip().upper() == 'U':
                edpos.append((affected_gene, mf.location, True))
    return edpos


def regexp(expr, item):
    """Regexp for search in sqlite database"""
    reg = re.compile(expr, re.IGNORECASE)
    if reg.search(item) is not None:
        return True
    else:
        exprl = [re.compile(ex, re.IGNORECASE) for ex in expr.split()]
        hasmatch = [x.search(item) for x in exprl]
        if all(hasmatch):
            return True
        else:
            # if written as x REGEXP y, item is x while expr is y
            # thereofre there is a need to escape item since it come from the
            # now use item as searching param
            iteml = [re.compile(re.escape(
                repr(itm)[1:-1]).replace('\\\\', '\\'), re.IGNORECASE) for itm in item.split()]
            hasmatch = [x.search(expr) for x in iteml]
            if all(hasmatch):
                return True
    return False


def collect_edpos_from_db(genome, gene, conn=None, dbfile=None, complete=False):
    """Collect known editing positions from a database directly"""
    if not dbfile and not conn:
        raise ValueError("You either need a connector or the database file")
    if not conn:
        conn = sqlite3.connect(dbfile)
        conn.create_function("REGEXP", 2, regexp)

    query = """SELECT accession, organism, location, gene, gstatus, evidence, type, 
                    nuc1, nuc2, pos, size, codchange, gcode, dnaseq, rnaseq, protseq
                    FROM editing WHERE (organism REGEXP ? or genbank REGEXP ?) and gene REGEXP ?"""
    if complete:
        query += " and gstatus='complete'"
    cur = conn.execute(query, [genome, genome, gene])
    rows = cur.fetchall()
    return rows


def get_relative_position(entry, ref_feat, sgene_reg, match_name=False, spec=None):
    """Return a gene specific position, from a genome location object.
    The offset of the relative position is 0"""
    gname, edloc, RNA_label = entry

    if (not match_name or not gname or sgene_reg.lower() == gname.lower()):
        dist = 0
        reflocs = ref_feat.location.parts
        strands = list(set([x.strand for x in reflocs]))
        #assert len(strands)==1, "Trans splicing in {} ({}): {}".format(sgene_reg, spec, strands)
        if not ref_feat.location.strand:
            MSG.add("==> Trans splicing in {} ({}). Consider removing gene !".format(
                sgene_reg, spec))

        strands = ref_feat.location.strand or strands[0]
        if strands < 0:
            reflocs = list(reversed(reflocs))
        for reflocpart in reflocs:
            start = reflocpart.start
            end = reflocpart.end
            dist += abs(end - start)
            if not (edloc.start < start or edloc.end > end):
                return dist - (end - edloc.end), RNA_label
    return None


@click.command(help="Build dataset of true positive for C-to-U RNA editing")
@click.argument('records', nargs=-1, type=click.Path(exists=True))
@click.option('--genelist', type=click.Path(exists=True), help='List of genes to keep in a file')
@click.option('--gcode', default=1, type=click.INT, help='Reference genetic code to use for all genomes.')
@click.option('--revgenes',  help='Map containing gene synonymes.')
@click.option('--nmatch',  is_flag=True, help='Whether RNA editing misc_feature should match the given gene name or not.')
@click.option('--dbcheck',  help='Check and extend list of edited positions with database reference (union not intersection)')
@click.option('--save_prot', is_flag=True, help='Save the true protein sequences after editing is applied')
@click.option('-v', '--verbose',  is_flag=True, help='Enable verbose')
@click.option('-nd', '--no_dup',  is_flag=True, help='Do not allow gene duplicates (e.g: cox1_2, cox1_1)')
@click.option('--ignore_syn',  is_flag=True, help='Ignore positions where editing result in synonymous codons')
@click.option('--permissive',  default=0, type=click.IntRange(0, 2, clamp=True), help='Try to use alignment consensus to guess missing RNA editing position. This is only performed at positions where at least one RNA editing is already reported. If set to 1, only allow RNA editing in genome that already have it at other position of the same gene !')
@click.option('-wd', '--outdir',  default="./dataset/", help='Output directory where to save the dataset and all results')
def cli(records, genelist=None, gcode=1, revgenes=None, nmatch=False, dbcheck=None, save_prot=False, verbose=False, no_dup=False, ignore_syn=False, permissive=0, outdir=None):
    """Extract genome content based on a list of species """
    spec_code_map = {}
    gene2spec = ddict(dict)
    prot2spec = ddict(dict)
    trans2spec = ddict(dict)
    knownedited = ddict(dict)

    def motif_in_seq(motif, slist, product):
        return sum([1 for x in slist if (motif in x.lower())]) > 0 and \
            'hypothetical' in " ".join(product).lower()

    # get the and decode the file of each genomes:
    revgenes = load_rev_genes(revgenes) if revgenes else {}
    gfile2rec = load_genbank_file(records)
    genelist = parse_genelist(genelist)
    dbconn = None
    if dbcheck:
        dbconn = sqlite3.connect(dbcheck)
        dbconn.create_function("REGEXP", 2, regexp)

    # fetch gene and protein from the database and add them to the list
    for (spec, cur_seq) in gfile2rec.items():
        g_found_in_spc = {}
        sname = cur_seq.name  # id of the sequence
        description = "_".join(cur_seq.annotations.get(
            'organism', cur_seq.description).split()).replace('.', '')

        putative_ed_pos = collect_edpos_from_seqrec(cur_seq)
        for pos, f in enumerate(cur_seq.features):
            # this should filter hypothetical protein that we do not want
            sgene = None
            meet_condition = f.type.upper() == 'CDS' and 'gene' in f.qualifiers.keys(
            ) and not motif_in_seq('orf', f.qualifiers['gene'], f.qualifiers.get('product', []))
            if meet_condition:
                sgene = f.qualifiers['gene'][0].lower()
                if not g_found_in_spc.get(sgene.split('-')[0]):
                    sgene = sgene.split('-')[0]

                sgene = revgenes.get(sgene, sgene).lower()
                if no_dup:
                    sgene = sgene.split('_')[0]
                if genematch(sgene, genelist):
                    matching_pos = [get_relative_position(
                        ed_entry, f, sgene, match_name=nmatch, spec=spec) for ed_entry in putative_ed_pos]
                    # excision to remove none values
                    # print(matching_pos)
                    seq = None
                    seq = f.extract(cur_seq.seq)
                    if len(seq) % 3 != 0:
                        # print("%s | %d"%(spec, len(seq)%3))
                        try:
                            polyaterm = f.qualifiers['transl_except'][0]
                            pos_range, aa_term = polyaterm.strip(
                                ')').strip('(').split(',')
                            pos_range = pos_range.strip().split(':')[-1]
                            aa_term = aa_term.strip().split(':')[-1]
                            # adding polyA to complete mRNA
                            if 'TERM' in aa_term.upper():
                                # Get the number of A to add
                                n_A = (len(pos_range.split('..')) % 2) + 1
                                seq = seq + Seq('A' * n_A, seq.alphabet)
                                assert len(seq) % 3 == 0
                                MSG.add("FIXED : partial termination for the following gene: %s - %s | %d %d A added" %
                                        (sgene, spec, len(seq), n_A))

                        except:
                            MSG.add("Possible frame-shifting in the following gene : %s - %s | %s ==> %d (%d)" %
                                    (sgene, spec, sname, len(seq), len(seq) % 3))
                    if f.strand and f.strand < 0:  # fix reverse sequence positioning for edited site
                        matching_pos = sorted(
                            [abs(len(seq)-x[0]) for x in matching_pos if x is not None])
                    else:
                        matching_pos = sorted(
                            [x[0]-1 for x in matching_pos if x])

                    matching_pos = list(matching_pos)

                    if 'N' in seq:
                        MSG.add("Sequence with undefined nucleotide : %s - %s | %d" %
                                (sgene, spec, len(seq)))
                    try:
                        table = int(f.qualifiers['transl_table'][0])
                        spec_code_map[sname] = table
                    except:
                        pass
                    rec = SeqRecord(seq, id=sname, name=sname,
                                    description=description)
                    transseq = seq.translate(
                        table=spec_code_map.get(sname, gcode))
                    protseq = Seq(f.qualifiers.get(
                        "translation", [transseq._data])[0])
                    if transseq.endswith('*') or len(protseq) * 3 == len(rec) - 3:
                        protseq._data = protseq._data.rstrip('*') + '*'

                    if dbconn:
                        rows = collect_edpos_from_db(
                            description, sgene, conn=dbconn, complete=True)
                        potential_extension = set([])
                        MSG.add("Match for %s: %s" %
                                (description, str(len(rows) > 0)))
                        # SELECT accession, organism, location, gene, gstatus, evidence, type,
                        #nuc1, nuc2, pos, size, codchange, gcode, dnaseq, rnaseq, protseq
                        # FROM editing"""
                        for row in rows:
                            curdnaseq = row[13]
                            if "experiment" in row[5].lower() and row[6] == "substitution" and curdnaseq == seq._data:
                                pos = row[9]
                                # db  sequences are indexed from 1
                                potential_extension.add(pos-1)

                        matching_pos.extend(potential_extension)
                        matching_pos = list(set(matching_pos))

                    # force changes into protseq
                    if protseq._data == transseq._data and len(matching_pos) > 0:
                        #  none of the editing positions was considered not considered
                        raw_prot = list(protseq._data)
                        for nucpos in matching_pos:
                            cod_pos = nucpos // 3
                            aa_at_pos = str(
                                rec[cod_pos*3:cod_pos*3 + 3].translate(table=spec_code_map.get(sname, gcode)).seq)
                            #print(cod_pos, nucpos, rec[cod_pos*3:cod_pos*3 + 3], aa_at_pos)
                            if aa_at_pos != raw_prot[cod_pos]:
                                raw_prot[cod_pos] = aa_at_pos
                        protseq._data = "".join(raw_prot)

                    # convert seq to recseq now
                    transseq = SeqRecord(
                        transseq, id=sname, name=sname, description=description)

                    protseq = SeqRecord(
                        protseq, id=sname, name=sname, description=description)

                    if g_found_in_spc.get(sgene, 0) < 1:
                        # this is to ensure that the same gene is not added
                        # multiple time
                        assert len(rec) == len(protseq)*3 == len(transseq)*3, '%s: %s %d, %d, %d Length of gene not matching: \n%s' % (
                            sname, sgene, len(protseq), len(transseq), len(rec), rec.format('fasta'))
                        gene2spec[sgene][sname] = rec
                        prot2spec[sgene][sname] = protseq
                        trans2spec[sgene][sname] = transseq
                        knownedited[sgene][sname] = matching_pos

                    g_found_in_spc[sgene] = 1
    #triple_format(gene2spec, prot2spec, trans2spec, glist=['cox1'])
    msa_data = {}
    for gname, slist in prot2spec.items():
        true_seq_aln = None
        if len(slist) > 1:
            msa_data[gname] = align(gname, slist.values(
            ), 'mafft', outdir=os.path.join(outdir, "prot_align/"))
            edited_pos, cod_align = codon_aligner(
                gname, msa_data[gname], gene2spec[gname], trans2spec[gname], knownedited[gname], spec_code_map, ignore_syn=ignore_syn)
            true_seq_aln = get_correct_seq(
                cod_align, edited_pos, gcode=spec_code_map)
            if permissive > 0:
                edited_pos = complete_from_consensus(
                    true_seq_aln, cod_align, edited_pos, gcode=spec_code_map, only_ed_G=(permissive == 1))
                true_seq_aln = get_correct_seq(
                    cod_align, edited_pos, gcode=spec_code_map)
            if not save_prot:
                true_seq_aln = None
            save_data(gname, cod_align.values(), edited_pos, os.path.join(
                outdir, "editing_map"), trueseq=true_seq_aln)
    if dbconn:
        dbconn.close()
    if verbose and MSG:
        for msg in MSG:
            logging.info(msg)
    return gene2spec, prot2spec, trans2spec, spec_code_map


if __name__ == '__main__':
    cli()
