import numpy as np
from collections import Counter
from Bio.Align import AlignInfo, MultipleSeqAlignment


def get_cost(costfn, res1, res2, defcost=0, ignore_gap=False, gap_char='-', gap_cost=None):
    """Return the cost between residue1 and residue2
    costfn can be a function that accepts 3 arguments. For example, if you want the edit distance:
    costfn = lambda x,y,z: (x!=y)*1
    Alternatively costfn can also be a dict object that specify the substitution cost between a pair of residues
    set ignore_gap if you want to ignore gap to gap alignment in the score (score set to )
    """
    res1, res2 = res1.upper(), res2.upper()
    score = None
    if ignore_gap and res1 == res2 == gap_char:
        score = defcost
    elif gap_char and res1 != res2 and gap_char in [res1, res2]:
        score = gap_cost
    if score is None:
        if callable(costfn):
            score = costfn(res1, res2, defcost)
        else:
            score = costfn.get(res1, {}).get(res2) or costfn.get(res2, {}).get(res1)
        score = defcost if score is None else score
    return score


def compute_SP_score(alignment, cols=None, costfn={}, count_method=None, defcost=0, ignore_gap=True, gap_char='-', gap_cost=None):
    """Compute SP score for a given col of the alignment"""
    if not isinstance(alignment, np.ndarray):
        alignment = np.asarray(alignment)
    if cols is None:
        cols = np.arange(alignment.shape[-1])
    elif not isinstance(cols, list):
        cols = [cols]
    sp_vec = np.zeros(len(cols))
    for icol, pos in enumerate(cols):
        sp = 0
        col = alignment[:, cols[pos]]
        if count_method:
            col = count_method(col)
        print(col)
        col = list(col)
        c_col = list(Counter(col).items())
        n_res = len(c_col)
        for i, (res, oc) in enumerate(c_col):
            sp += (oc - 1)*get_cost(costfn, res, res,
                   defcost, ignore_gap, gap_char, gap_cost)
            for j in range(i, n_res):
                res2, oc2 = c_col[j]
                sp += oc*oc2*get_cost(costfn, res, res2, defcost,
                                      ignore_gap, gap_char, gap_cost)
        sp_vec[icol] = sp
    return sp_vec


def compute_ic_content(alignment, cols):
    """Compute the information content for the list of cols specified"""
    msa = alignment
    if isinstance(msa, np.ndarray):
        msa = MultipleSeqAlignment(msa)
    if not isinstance(msa, MultipleSeqAlignment):
        raise ValueError("Expect MultipleSeqAlignment or np.ndarray, got {}".format(msa.__class__))
   
    align_info = AlignInfo.SummaryInfo(msa)
    align_info.information_content()
    ic_vector = align_info.ic_vector
    # biopython wrong version hack
    if isinstance(ic_vector, dict):
        ic_vector = np.zeros(len(align_info.ic_vector))
        for (ic_i, ic_v) in align_info.ic_vector.items():
            ic_vector[ic_i] = ic_v
    return list(ic_vector)


def information_content(self, start=0,
                        end=None,
                        e_freq_table=None, log_base=2,
                        chars_to_ignore=None, pseudo_count=0):
    """Calculate the information content for each residue along an alignment.
    Arguments:
     - start, end - The starting an ending points to calculate the
       information content. These points should be relative to the first
       sequence in the alignment, starting at zero (ie. even if the 'real'
       first position in the seq is 203 in the initial sequence, for
       the info content, we need to use zero). This defaults to the entire
       length of the first sequence.
     - e_freq_table - A FreqTable object specifying the expected frequencies
       for each letter in the alphabet we are using (e.g. {'G' : 0.4,
       'C' : 0.4, 'T' : 0.1, 'A' : 0.1}). Gap characters should not be
       included, since these should not have expected frequencies.
     - log_base - The base of the logathrim to use in calculating the
       information content. This defaults to 2 so the info is in bits.
     - chars_to_ignore - A listing of characters which should be ignored
       in calculating the info content. Defaults to none.
    Returns:
     - A number representing the info content for the specified region.
    Please see the Biopython manual for more information on how information
    content is calculated.
    """
    # if no end was specified, then we default to the end of the sequence
    if end is None:
        end = len(self.alignment[0].seq)
    if chars_to_ignore is None:
        chars_to_ignore = []

    if start < 0 or end > len(self.alignment[0].seq):
        raise ValueError("Start (%s) and end (%s) are not in the \
                range %s to %s"
                         % (start, end, 0, len(self.alignment[0].seq)))
    # determine random expected frequencies, if necessary
    random_expected = None
    if not e_freq_table:
        # TODO - What about ambiguous alphabets?
        base_alpha = Alphabet._get_base_alphabet(self.alignment._alphabet)
        if isinstance(base_alpha, Alphabet.ProteinAlphabet):
            random_expected = Protein20Random
        elif isinstance(base_alpha, Alphabet.NucleotideAlphabet):
            random_expected = Nucleotide4Random
        else:
            errstr = "Error in alphabet: not Nucleotide or Protein, "
            errstr += "supply expected frequencies"
            raise ValueError(errstr)
        del base_alpha
    elif not isinstance(e_freq_table, FreqTable.FreqTable):
        raise ValueError("e_freq_table should be a FreqTable object")

    # determine all of the letters we have to deal with
    all_letters = self._get_all_letters()
    for char in chars_to_ignore:
        all_letters = all_letters.replace(char, '')

    info_content = {}
    for residue_num in range(start, end):
        freq_dict = self._get_letter_freqs(residue_num,
                                           self.alignment,
                                           all_letters,
                                           chars_to_ignore,
                                           pseudo_count,
                                           e_freq_table,
                                           random_expected)
        # print freq_dict,
        column_score = self._get_column_info_content(freq_dict,
                                                     e_freq_table,
                                                     log_base,
                                                     random_expected)
        info_content[residue_num] = column_score
    # sum up the score
    total_info = sum(info_content.values())
    # fill in the ic_vector member: holds IC for each column
    # reset ic_vector to empty list at each call
    self.ic_vector = []
    for (i, k) in enumerate(info_content):
        self.ic_vector.append(info_content[i + start])
    return total_info

def _get_letter_freqs(self, residue_num, all_records, letters, to_ignore,
                      pseudo_count=0, e_freq_table=None, random_expected=None):
    """Determine the frequency of specific letters in the alignment (PRIVATE).
    Arguments:
     - residue_num - The number of the column we are getting frequencies
       from.
     - all_records - All of the SeqRecords in the alignment.
     - letters - The letters we are interested in getting the frequency
       for.
     - to_ignore - Letters we are specifically supposed to ignore.
     - pseudo_count - Optional argument specifying the Pseudo count (k)
       to add in order to prevent a frequency of 0 for a letter.
     - e_freq_table - An optional argument specifying the expected
       frequencies for each letter. This is a SubsMat.FreqTable instance.
     - random_expected - Optional argument that specify the frequency to use
       when e_freq_table is not defined.
    This will calculate the frequencies of each of the specified letters
    in the alignment at the given frequency, and return this as a
    dictionary where the keys are the letters and the values are the
    frequencies. Pseudo count can be added to prevent a null frequency
    """
    freq_info = self._get_base_letters(letters)

    total_count = 0

    gap_char = self._get_gap_char()

    if pseudo_count < 0:
        raise ValueError("Positive value required for "
                         "pseudo_count, %s provided" % (pseudo_count))

    # collect the count info into the dictionary for all the records
    for record in all_records:
        try:
            if record.seq[residue_num] not in to_ignore:
                weight = record.annotations.get('weight', 1.0)
                freq_info[record.seq[residue_num]] += weight
                total_count += weight
        # getting a key error means we've got a problem with the alphabet
        except KeyError:
            raise ValueError("Residue %s not found in alphabet %s"
                             % (record.seq[residue_num],
                                self.alignment._alphabet))

    if e_freq_table:
        if not isinstance(e_freq_table, FreqTable.FreqTable):
            raise ValueError("e_freq_table should be a FreqTable object")

        # check if all the residus in freq_info are in e_freq_table
        for key in freq_info:
            if (key != gap_char and key not in e_freq_table):
                raise ValueError("letters in current column %s "
                                 "and not in expected frequency table %s"
                                 % (list(freq_info) - [gap_char],
                                    list(e_freq_table)))

    if total_count == 0:
        # This column must be entirely ignored characters
        for letter in freq_info:
            assert freq_info[letter] == 0
            # TODO - Map this to NA or NaN?
    else:
        # now convert the counts into frequencies
        for letter in freq_info:
            if pseudo_count and (random_expected or e_freq_table):
                # use either the expected random freq or the
                if e_freq_table:
                    ajust_freq = e_freq_table[letter]
                else:
                    ajust_freq = random_expected

                ajusted_letter_count = freq_info[
                    letter] + ajust_freq * pseudo_count
                ajusted_total = total_count + pseudo_count
                freq_info[letter] = ajusted_letter_count / ajusted_total

            else:
                freq_info[letter] = freq_info[letter] / total_count

    return freq_info

def _get_column_info_content(self, obs_freq, e_freq_table, log_base,
                             random_expected):
    """Calculate the information content for a column (PRIVATE).
    Arguments:
     - obs_freq - The frequencies observed for each letter in the column.
     - e_freq_table - An optional argument specifying the expected
       frequencies for each letter. This is a SubsMat.FreqTable instance.
     - log_base - The base of the logathrim to use in calculating the
       info content.
    """
    gap_char = self._get_gap_char()

    if e_freq_table:
        if not isinstance(e_freq_table, FreqTable.FreqTable):
            raise ValueError("e_freq_table should be a FreqTable object")
        # check the expected freq information to make sure it is good
        for key in obs_freq:
            if (key != gap_char and key not in e_freq_table):
                raise ValueError("Expected frequency letters %s "
                                 "do not match observed %s"
                                 % (list(e_freq_table),
                                    list(obs_freq) - [gap_char]))

    total_info = 0.0

    for letter in obs_freq:
        inner_log = 0.0
        # if we have expected frequencies, modify the log value by them
        # gap characters do not have expected frequencies, so they
        # should just be the observed frequency.
        if letter != gap_char:
            if e_freq_table:
                inner_log = obs_freq[letter] / e_freq_table[letter]
            else:
                inner_log = obs_freq[letter] / random_expected
        # if the observed frequency is zero, we don't add any info to the
        # total information content
        if inner_log > 0:
            letter_info = (obs_freq[letter] *
                           math.log(inner_log) / math.log(log_base))
            total_info += letter_info
    return total_info


def get_column(self, col):
    # TODO - Deprecate this and implement slicing?
    return self.alignment[:, col]
