# More than 80% of the code here was copied from coretracker/coreutils/Faces.py
# which was written for displaying a tree with a list of sequences

import argparse
import math

from PyQt5.QtCore import Qt, QPointF
from PyQt5.QtWidgets import (
    QGraphicsRectItem, QGraphicsLineItem, QGraphicsSimpleTextItem, QGraphicsTextItem)
from PyQt5.QtGui import (QPen, QLinearGradient,
                         QColor, QBrush, QFont)
from ete3 import faces, Tree, TreeStyle, NodeStyle
from ete3.treeview.main import COLOR_SCHEMES
from collections import OrderedDict, defaultdict as ddict
import os
import sys
from Bio import AlignIO, SeqIO
from Bio.Data import CodonTable
import logging
from ete3 import TextFace, AttrFace

# change this to another color if you want to set how node with at least one mark are displayed
MARKED_NODE_COLOR = "red"
RES_COL_WIDTH = 14

_aafgcolors = {
    'A': "#212121",
    'R': "#212121",
    'N': "#212121",
    'D': "#212121",
    'C': "#212121",
    'Q': "#212121",
    'E': "#212121",
    'G': "#212121",
    'H': "#212121",
    'I': "#212121",
    'L': "#212121",
    'K': "#212121",
    'M': "#212121",
    'F': "#212121",
    'P': "#212121",
    'S': "#212121",
    'T': "#212121",
    'W': "#212121",
    'Y': "#212121",
    'V': "#212121",
    'B': "#212121",
    'Z': "#212121",
    'X': "#212121",
    '.': "#212121",
    '-': "#212121",
    '*': "#212121",
}

_aabgcolors = {
    'A': "#C8C8C8",
    'R': "#145AFF",
    'N': "#00DCDC",
    'D': "#E60A0A",
    'C': "#E6E600",
    'Q': "#00DCDC",
    'E': "#E60A0A",
    'G': "#DBDBDB",
    'H': "#8282D2",
    'I': "#0F820F",
    'L': "#0F820F",
    'K': "#145AFF",
    'M': "#E6E600",
    'F': "#4e4eba",
    'P': "#DC9682",
    'S': "#FA9600",
    'T': "#FA9600",
    'W': "#B45AB4",
    'Y': "#4e4eba",
    'V': "#0F820F",
    'B': "#FF69B4",
    'Z': "#FF69B4",
    'X': "#BEA06E",
    '.': "#FFFFFF",
    '-': "#FFFFFF",
    '*': "#FFFFFF",
}

_ntfgcolors = {
    'A': '#000000',
    'G': '#000000',
    'I': '#000000',
    'C': '#000000',
    'T': '#000000',
    'U': '#000000',
    '.': "#000000",
    '-': "#000000",
    ' ': "#000000"
}

_ntbgcolors = {
    'A': '#A0A0FF',
    'G': '#FF7070',
    'I': '#80FFFF',
    'C': '#FF8C4B',
    'T': '#A0FFA0',
    'U': '#FF8080',
    '.': "#FFFFFF",
    '-': "#FFFFFF",
    ' ': "#FFFFFF"
}


def _get_codon_fgcolors(codontable):
    """Get colon foreground color"""
    return dict((k, _aafgcolors[v]) for (k, v) in codontable.items())


def _get_codon_bgcolors(codontable, spec_codon_col):
    """Get colon background color"""
    return dict((k, spec_codon_col.get(k, _aabgcolors.get(v, "#FFFFFF"))) for (k, v) in codontable.items())


class _LineItem(QGraphicsLineItem):

    def __init__(self, w, h, fgcolor):
        QGraphicsLineItem.__init__(self)
        self.setLine(w / 2., 0, w / 2., h)
        if fgcolor:
            self.setPen(QPen(QColor(fgcolor)))
        else:
            self.setPen(QPen(QColor('#000000')))


class LineFace(faces.Face):
    """
    Creates a Line face.
    """

    def __init__(self, width, height, fgcolor):
        faces.Face.__init__(self)
        self.width = width
        self.height = height
        self.fgcolor = fgcolor
        self.type = "item"
        self.rotable = True

    def update_items(self):
        self.item = _LineItem(self.width, self.height, self.fgcolor)

    def _width(self):
        return self.width

    def _height(self):
        return self.height


class _LineItem(QGraphicsLineItem):

    def __init__(self, w, h, fgcolor):
        QGraphicsLineItem.__init__(self)
        self.setLine(w / 2., 0, w / 2., h)
        if fgcolor:
            self.setPen(QPen(QColor(fgcolor)))
        else:
            self.setPen(QPen(QColor('#000000')))


class SequenceFace(faces.StaticItemFace):
    """
    Creates a new molecular sequence face object.
    :param seq: Sequence string to be drawn
    :param seqtype: Type of sequence: "nt" or "aa"
    :param fsize: Font size, (default=10)
    You can set custom colors for amino-acids or nucleotides:
    :param None codon: a string that corresponds to the reverse
      translation of the amino-acid sequence
    :param None col_w: width of the column (if col_w is lower than
      font size, letter wont be displayed)
    :param None fg_colors: dictionary of colors for foreground, with
      as keys each possible character in sequences, and as value the
      colors
    :param None bg_colors: dictionary of colors for background, with
      as keys each possible character in sequences, and as value the
      colors
    """

    def __init__(self, seq, seqtype="aa", fsize=10, ftype="Monospace",
                 fg_colors=None, bg_colors={}, col_w=None, codontable={}, black_out=[]):
        self.seq = seq
        self.fsize = fsize
        self.style = seqtype
        self.col_w = float(self.fsize + 0.5) if col_w is None else float(col_w)
        self.width = 0  # will store the width of the whole sequence
        black_out = sorted(black_out)
        self.mark = {}
        self.ftype = ftype
        multipliers = 1

        if self.style == "aa":
            if not fg_colors:
                fg_colors = _aafgcolors
            if not bg_colors:
                bg_colors = _aabgcolors

        elif self.style == 'codon':
            multipliers = 3
            self.col_w *= 3
            if not isinstance(self.seq, list):
                # only consider the position where 3 nuc can be obtained
                self.seq = [self.seq[i:i + 3] for i in range(0,
                                                             len(self.seq) - len(self.seq) % 3, 3)]
            if not fg_colors:
                fg_colors = _get_codon_fgcolors(codontable)

            bg_colors.update(_get_codon_bgcolors(codontable, bg_colors))

        else:
            if not fg_colors:
                fg_colors = _ntfgcolors
            if not bg_colors:
                bg_colors = _ntbgcolors

        def __init_colors(color_dic, keys, col='#212121', retcolor=False):
            """to speed up the drawing of colored rectangles and characters"""
            new_color_dic = {}
            for residue in keys:
                # use default black color if color not found
                c = color_dic.get(residue, col)
                if retcolor:
                    new_color_dic[residue] = c
                else:
                    new_color_dic[residue] = QBrush(QColor(c))
            return new_color_dic

        self.fg_col = __init_colors(fg_colors, set(self.seq), retcolor=True)
        self.bg_col = __init_colors(bg_colors, set(self.seq), '#FFFFFF')

        for pos, residue in enumerate(self.seq):
            min_pos, max_pos = pos*multipliers, (pos+1)*multipliers
            for mpos in black_out:
                if min_pos <= mpos and mpos < max_pos:
                    # save position in residue to mark
                    self.mark[pos] = mpos % multipliers

        # for future?
        self.row_h = 14.0
        super(SequenceFace, self).__init__(None)

    def update_items(self):
        rect_cls = QGraphicsRectItem
        self.item = rect_cls(0, 0, self.width, self.row_h)
        seq_width = 0
        nopen = QPen(Qt.NoPen)
        font = QFont(self.ftype, self.fsize)
        for i, letter in enumerate(self.seq):
            width = self.col_w
            rectitem = rect_cls(0, 0, width, self.row_h, parent=self.item)
            rectitem.setX(seq_width)  # to give correct X to children item
            rectitem.setBrush(self.bg_col[letter])
            rectitem.setPen(nopen)
            # write letter if enough space
            if width >= self.fsize:
                text = None
                if self.style == "codon":
                    text = QGraphicsTextItem(parent=rectitem)
                    text.setHtml(self.get_html(i, self.fg_col[letter]))
                else:
                    text = QGraphicsSimpleTextItem(letter, parent=rectitem)
                    text.setBrush(QBrush(QColor(self.fg_col[letter])))
                text.setFont(font)
                # Center text according to rectitem size
                txtw = text.boundingRect().width()
                txth = text.boundingRect().height()
                text.setPos((width - txtw) / 2, (self.row_h - txth) / 2)
            seq_width += width
        self.width = seq_width

    def get_html(self, pos, fgcolor):
        letter = self.seq[pos]
        if not self.mark.get(pos):
            return '<span style="color:%s;">%s</span>'%(fgcolor, letter)
        else:
            html = ""
            for i, l in enumerate(letter):
                if i != self.mark[pos]:
                    html += l
                else:
                    html += '<span style="color:%s; background-color:%s;margin-left:1px; margin-right:1px;">%s</span>' % (
                        "#FFFFFF", "#111111", l)
            return '<span style="color:%s;">%s</span>'%(fgcolor, html)


class List90Face(faces.StaticItemFace):
    """Static text Face object
    :param l:        List of element to be drawn
    :param fsize:    Font size, e.g. 10,12,6, (default=10)
    :param fgcolor:  Foreground font color. RGB code or color name in :data:`SVG_COLORS`
    :param bgcolor:  Background font color. RGB code or color name in :data:`SVG_COLORS`
    """

    def __init__(self, l, ftype="Courier", fstyle="normal", fsize=10,
                 fgcolor="black", bgcolor="white", col_w=14.0, rotation=90):
        self.liste = l
        self.ftype = ftype
        self.fgcolor = fgcolor
        self.bgcolor = bgcolor
        self.fsize = fsize
        self.row_h = float(self.fsize + 1)
        self.col_w = col_w
        self.width = 0
        self.rot = rotation
        self.fstyle = fstyle
        self.coeff_h = max([len(str(x)) for x in self.liste])

        super(List90Face, self).__init__(None)

    def __repr__(self):
        return "Text Face [%s] (%s)" % (self._text, hex(self.__hash__()))

    def get_text(self):
        return self._text

    def update_items(self):
        self.item = QGraphicsRectItem(
            0, 0, self.width, self.row_h * self.coeff_h)
        seq_width = 0
        nopen = QPen(Qt.NoPen)
        self.item.setPen(nopen)
        font = QFont(self.ftype, self.fsize)
        if self.fstyle == "italic":
            font.setStyle(QFont.StyleItalic)
        elif self.fstyle == "oblique":
            font.setStyle(QFont.StyleOblique)
        rect_cls = QGraphicsRectItem
        for i, val in enumerate(self.liste):
            width = self.col_w
            height = self.row_h * len(str(val)) + 1
            rectitem = rect_cls(0, 0, width, height, parent=self.item)
            rectitem.setX(seq_width)  # to give correct X to children item
            rectitem.setBrush(QBrush(QColor(self.bgcolor)))
            rectitem.setPen(nopen)

            # write letter if enough space in height
            if height >= self.fsize:
                text = QGraphicsSimpleTextItem(str(val), parent=rectitem)
                text.setFont(font)
                text.setBrush(QBrush(QColor(self.fgcolor)))
                # Center text according to rectitem size
                txtw = text.boundingRect().width() /3.0
                txth = text.boundingRect().height()
                text.setRotation(self.rot)
                text.setX(txth*1.5)
                #text.setY(0)
            seq_width += width
        self.width = seq_width


def format_tree(tree, alignment, al_len_dict, edpos, codontable={}, colors=None, codon_col={}, text="C-to-U RNA editing", ic_contents=[]):
    """Format the rendering of tree data for alignment"""
    t = tree.copy()
    # alignment is ordered dict

    # flip alignment dict from gene ==> species ==> seq
    # to species ==> gene ==> seq
    specSeq = ddict(str)
    edposSeq =  ddict(list)
    cur_len = 0
    limits = []
    for gname, specdict in alignment.items():
        for node in t:
            # fill missing with gap
            specSeq[node.name] += specdict.get(node.name,
                                               al_len_dict[gname]*'-')
            edposSeq[node.name] += [
                x+cur_len for x in edpos[gname].get(node.name, [])]
            # if node.name == 'Y08501':
            #     print(gname)
            #     print( edposSeq[node.name])
        cur_len += al_len_dict.get(gname, 0)
        limits.append((gname, cur_len))

    for node in t:
        node.add_feature("sequence", specSeq[node.name])
        node.add_feature('edlist', edposSeq[node.name])

    ts = TreeStyle()
    ts.branch_vertical_margin = 15
    ts.scale = 15
    ts.allow_face_overlap = False
    ts.show_scale = False
    ts.show_leaf_name = False

    ns = NodeStyle()
    ns['shape'] = 'square'
    ns['fgcolor'] = 'black'
    ns['size'] = 0

    def layout(node):
        node.img_style = ns
        if node.is_leaf():
            faces.add_face_to_node(
                AttrFace('fullname', fsize=14, fgcolor=(MARKED_NODE_COLOR if (node.name in colors or node.fullname in colors) else 'black')), node, 0, position="aligned")
            if hasattr(node, "sequence") and node.sequence:
                seqface = SequenceFace(
                    node.sequence, "codon", fsize=13, codontable=codontable,  col_w=RES_COL_WIDTH,  bg_colors=codon_col, black_out=node.edlist)
                faces.add_face_to_node(seqface, node, 1, position="aligned")

    ts.layout_fn = layout

    # ts.title.add_face(TextFace('(%s) - SP score : %.0f | IC = %.2f' % (codon, sum(SP_score), sum(ic_contents)),
    #                            fsize=14, fgcolor='red'), 0)
    # ts.aligned_header.add_face(
    #     faces.RectFace(14, 14, 'white', 'white'), 1)

    # ts.aligned_foot.add_face(
    #     faces.RectFace(14, 14, 'white', 'white'), 1)

    # for (cod, col) in codon_col.items():
    #     ts.legend.add_face(faces.RectFace(50, 25, col, col), column=0)
    #     ts.legend.add_face(TextFace("  %s " % cod, fsize=8), column=1)

    ts.legend.add_face(TextFace(text, fsize=14), column=1)
    ts.legend_position = 1

    ind = 1
    prev_gend = 0
    for (gname, gend) in limits:
        ts.aligned_foot.add_face(List90Face(list(range(0, gend - prev_gend, 3)), fsize=13, ftype='Monospace', col_w=RES_COL_WIDTH*3), ind)
        ts.aligned_foot.add_face(faces.RectFace(
            RES_COL_WIDTH*(gend - prev_gend), 13, '#BBBBBB', '#EEEEEE'), ind)
        ts.aligned_foot.add_face(TextFace(gname, fsize=13), ind)
        ts.aligned_foot.add_face(faces.RectFace(
            RES_COL_WIDTH*(gend - prev_gend), 5, 'white', 'white'), ind)
        prev_gend += gend
        ind += 1

    #t.dist = 0
    ts.margin_left = 5
    ts.margin_right = 5
    ts.margin_bottom = 5
    return t, ts


def display_tree(tree, aligndict, al_len_dict, eddict, outfile, gcode, colormap, dpi=1200, width=500):

    # check if leave name and seq name match
    gcode['---'] = '-'
    colors = set([y for x in eddict.values() for y in x.keys() if x[y]])
    t, ts = format_tree(tree, aligndict, al_len_dict, eddict, codontable=gcode,
                        codon_col=colormap, colors=colors, text="C-to-U RNA editing")
    #ts.title.add_face(TextFace("Alignment display ", fsize=14), column=0)
    t.render(outfile, dpi=dpi, tree_style=ts, w=width)

def get_ed_pos(edpath):
    ed = {}
    try:
        ED = open(edpath)
        for line in ED:
            line = line.strip()
            if line:
                line = line.split(':')
                ed[line[0].strip()] = [int(x.strip())
                                       for x in line[-1].strip().split(",")]
    except OSError:
        pass
    return ed


if __name__ == '__main__':
    im_choices = ("svg", "pdf", "png", "html")
    parser = argparse.ArgumentParser(
        "EdAlignViz", description='Display a codon-by-codon alignment with a phylogenetic tree, and show the list of edited position')
    parser.add_argument('-t', '--tree', required=True,
                        dest='tree', help="Path to the phylogenetic tree")
    parser.add_argument('-a', '--align', nargs="+", required=True,
                        dest='align', help="List of sequence alignments")
    parser.add_argument('-m', '--mark_ext', nargs="?", const="truth", dest='mark_ext',
                        help="The extension of the file in which position to mark are saved. If it is provided without argument, a '.truth' extension is used. Note that the basename  of the file should match the one of the respective alignment.")
    parser.add_argument('-o', '--output', default="output", dest='output',
                        help="Path to the output file where to save the plot")
    parser.add_argument('--fmt',  choices=im_choices,
                        default="svg", help="Image output format")
    parser.add_argument('--dpi',  type=int, default=300, help="DPI for image quality")
    parser.add_argument('-w','--width',  type=float, default=10000, help="Width of the image")

    parser.add_argument('--gcode', default=1, type=int,
                        help="Genetic code table to use. Only accept Standard table (NCBI) as for now")
    parser.add_argument('-v', '--verbose', dest="verbose",
                        action="store_true", help="Whether debugging text should be printed")
    parser.add_argument(
        '--colormap', help="Name of a file that map a codon to a specific hex color used as background. Space separated values")

    args = parser.parse_args()
    colormap = {}
    if args.colormap:
        with open(args.colormap) as COLOR:
            for line in COLOR:
                line = line.strip().split()
                if line:
                    colormap[line[0]] = line[1]

    out = args.output.rsplit(".", 1)
    if len(out) == 1:
        out = out[0] + "." + args.fmt
    elif len(out) == 2 and out[-1].lower() not in im_choices:
        out = out[0] + ".svg"
    else:
        out = out[0] + "."+out[-1].lower()

    tree = Tree(args.tree)
    al_per_gene = {}
    eddict = {}
    al_len = {}
    gcode = CodonTable.unambiguous_dna_by_id[args.gcode]
    gcode = gcode.forward_table
    spec_and_descr = {}
    for al in args.align:
        if os.path.isfile(al):
            try:
                gname, ext = os.path.splitext(al)
                gene = os.path.basename(gname)
                cur_msa = AlignIO.read(al, "fasta")
                # add alignment length
                al_len[gene] = cur_msa.get_alignment_length()
                spec_and_descr.update(dict((x.id, x.description)
                                           for x in cur_msa))
                # add sequences for each gene
                al_per_gene[gene] = dict((x.id, str(x.seq)) for x in cur_msa)
                # get edited positions in each gene
                if args.mark_ext:
                    ed = get_ed_pos(gname+"."+args.mark_ext.lstrip(","))
                    if ed:
                        eddict[gene] = ed
            except Exception as e:
                print(e)
                logging.warning(
                    "File %s cannot be read as fasta alignment. Skipping the file" % al)

    # now try to reset the leaves name in the tree
    common_name = []
    for leaf in tree:
        new_id = [x for x, y in spec_and_descr.items() if leaf.name in y]
        if new_id:
            common_name.append(new_id[0])
            leaf.add_feature("fullname", leaf.name)
            leaf.name = new_id[0]
    tree.prune(common_name)
    #print(tree.get_ascii(attributes=["name", "fullname"]))
    if len(tree) < len(spec_and_descr):
        logging.warning("The following genomes are missing in the tree:")

    for s in set.symmetric_difference(set(tree.get_leaf_names()), spec_and_descr.keys()):
        logging.warning(spec_and_descr[s])
    display_tree(tree, al_per_gene, al_len, eddict, out, gcode, colormap, dpi=args.dpi, width=args.width)