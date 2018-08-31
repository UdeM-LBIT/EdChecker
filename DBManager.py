#!/usr/bin/env python
import sqlite3
import os
import sys
import click
import itertools as it
import re

reSUB = re.compile("(.*)-->(.*)")  # match nucleotide conversion
reDEL = re.compile("<-(.*)->")  # match nucleotide deletion
reINS = re.compile("->(.*)<-")  # match nucleotide insertion


def create_db(filename):
    """Create a new database to store known RNA editing with Experimental evidence"""

    schema = '''
        create table IF NOT EXISTS editing (
            id integer primary key autoincrement not null,
            accession text not null,
            organism text not null,
            location text,
            genbank text,
            gtype text,
            gene text not null,
            gstatus text default 'complete',
            pubmed text,
            comment text,
            type text,
            evidence text,
            nuc1 text,
            nuc2 text,
            pos int not null,
            size int not null default 1,
            codchange text,
            dnaseq text,
            rnaseq text,
            protseq text,
            gcode int default 1
        );
        '''
    # this will be re-designed if needed
    with sqlite3.connect(filename) as conn:
        conn.executescript(schema)


def parse_data(lines):
    """Parse each entry of the raw editing database dump"""
    data = {}
    val = ""
    edlist = []
    sentinel = object()
    l = next(lines, sentinel)
    while l and (l is not sentinel):
        l = l.strip()
        if l.startswith('ACCESSION'):
            val = l.split()[1]
            data['accession'] = val
            l = next(lines, sentinel)
        elif l.startswith('ORGANISM'):
            val = l.split(None, 1)[1]
            data['organism'] = val
            l = next(lines, sentinel)
        elif l.startswith('LOCATION'):
            val = l.split()[1].lower()
            data['location'] = val
            l = next(lines, sentinel)
        elif l.startswith('SOURCE'):
            val = l.split()
            data['gtype'] = val[1]
            data['gene'] = val[2]
            if len(val) == 4:
                data['gstatus'] = val[3]
            l = next(lines, sentinel)
        elif l.startswith('GENBANK'):
            val = l.split(None, 1)
            if len(val) > 1:
                data["genbank"] = val[-1]
            l = next(lines, sentinel)
            while (l is not sentinel) and l.startswith(" "):
                data["genbank"] = data["genbank"] + " " + l.strip()
                l = next(lines, sentinel)
        elif l.startswith('PUBMED'):
            val = l.split(None, 1)
            if len(val) > 1:
                data["pubmed"] = val[-1]
            l = next(lines, sentinel)
            while (l is not sentinel) and l.startswith(" "):
                data["pubmed"] = data["pubmed"] + " " + l.strip()
                l = next(lines, sentinel)
        elif l.startswith('COMMENT'):
            val = l.split(None, 1)
            if len(val) > 1:
                data["comment"] = val[-1]
            l = next(lines, sentinel)
            while (l is not sentinel) and l.startswith(" "):
                data["comment"] = data["comment"] + " " + l.strip()
                l = next(lines, sentinel)
        elif l.startswith('EDITING'):
            val = l.split()
            l = next(lines, sentinel)
            position_reading = False
            eddata = {}
            eddata["type"] = val[-1].lower()
            edposdt = []
            while not (l is sentinel) and l.startswith(" ") and not l.startswith('ORIGIN'):
                l = l.strip()
                if l.startswith('DETAILS'):
                    val = l.split()
                    edtype = val[1]
                    evidence = ""
                    if len(val) > 2:
                        evidence = " ".join(val[2:])
                    eddata["evidence"] = evidence
                    if reSUB.match(edtype) or eddata.get('type').lower() == 'substitution':
                        eddata["nuc1"], eddata["nuc2"] = reSUB.match(
                            edtype).groups()
                    elif reINS.match(edtype) or eddata.get('type').lower() == 'insertion':
                        eddata["nuc1"] = reINS.match(edtype).groups()[0]
                    elif reDEL.match(edtype) or eddata.get('type').lower() == 'deletion':
                        eddata["nuc1"] = reDEL.match(edtype).groups()[0]
                if l.startswith("POSITIONS"):
                    position_reading = True

                elif position_reading and l:
                    l = l.split()
                    pos = l[0]
                    det = ""
                    if len(l) > 1:
                        det = " ".join(l[1:])
                    pos = int(pos.strip('*'))
                    size = 1
                    codchange = None
                    if eddata['type'] != 'substitution':
                        size = int(det.replace("EXT", "").strip())
                    else:
                        codchange = ":".join([x.strip()
                                              for x in det.split('-->')])
                    edposdt.append((pos, size, codchange))
                l = next(lines, sentinel)
            for ed in edposdt:
                edt = eddata.copy()
                edt["pos"] = ed[0]
                edt["size"] = ed[1]
                edt["codchange"] = ed[2]
                edlist.append(edt)

        elif l.startswith("ORIGIN"):
            l = next(lines, sentinel)
            while not (l is sentinel) and l.startswith(" "):
                if 'GENOMIC' in l.upper():
                    geneseq = ""
                    l = next(lines, sentinel)
                    while not (l is sentinel) and l.strip() not in ['GENOMIC', "cDNA", "PROTEIN"]:
                        if "BASE" not in l and "LENGTH" not in l:
                            geneseq += "".join(l.strip().split()[1:]).upper()
                        l = next(lines, sentinel)
                    data["dnaseq"] = geneseq

                elif 'CDNA' in l.upper():
                    geneseq = ""
                    l = next(lines, sentinel)
                    while not (l is sentinel) and l.strip() not in ['GENOMIC', "cDNA", "PROTEIN"]:
                        if "BASE" not in l and "LENGTH" not in l:
                            geneseq += "".join(l.strip().split()[1:]).upper()
                        l = next(lines, sentinel)
                    data["rnaseq"] = geneseq

                elif 'PROTEIN' in l.upper():
                    geneseq = ""
                    l = next(lines, sentinel)
                    while not (l is sentinel) and l.strip() not in ['GENOMIC', "cDNA", "PROTEIN"]:
                        if "TRANSLATION TABLE" not in l and "LENGTH" not in l:
                            geneseq += "".join(l.strip().split()[1:]).upper()
                        elif "TRANSLATION TABLE" in l:
                            data["gcode"] = int(l.strip().rsplit(None, 1)[-1])
                        l = next(lines, sentinel)
                    data["protseq"] = geneseq
    return data, edlist


@click.command(help="Create a new database that contains ground truth for codon identity")
@click.argument('infile', type=click.Path(exists=True))
@click.option('-dbout', default="./dataset/Edit.db", help='Outfile name for the database')
def db_fill(infile, dbout, verbose=False):
    """Create a database and save all editing entry inside the database"""
    create_db(dbout)
    i = 0
    with sqlite3.connect(dbout) as conn:
        query = "INSERT INTO editing ({}) VALUES ({})"
        with open(infile) as INDB:
            for key, line in it.groupby(INDB, key=lambda l: l.startswith('//Stop')):
                if not key:
                    i += 1
                    data, edlist = parse_data(line)
                    for ed in edlist:
                        dt = data.copy()
                        dt.update(ed)
                        # print(dt)
                        columns = ', '.join(dt.keys())
                        pholders = ', :'.join(dt.keys())
                        pholders = ':'+pholders
                        conn.execute(query.format(columns, pholders), dt)
                    # print("------------------------------------------")
                    if verbose:
                        print("Inserting entry: {}".format(i))

        conn.commit()
        print("{} entries read".format(i))
        print("{} elements inserted in {}".format(conn.total_changes, dbout))


if __name__ == '__main__':
    db_fill()
