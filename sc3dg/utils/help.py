help = {
    'select':'pairtools select CONDITION : A Python expression; if it returns True, select the read pair. \
             Any column declared in the #columns line of the pairs header can be accessed \
            by its name. If the header lacks the #columns line, the columns are assumed \
            to follow the .pairs/.pairsam standard (readID, chrom1, chrom2, pos1, pos2, \
            strand1, strand2, pair_type). Finally, CONDITION has access to COLS list\
            which contains the string values of columns. In Bash, quote CONDITION with\
            single quotes, and use double quotes for string variables inside CONDITION.',
    'add_columns':'Report extra columns describing alignments \
                    Possible values (can take multiple values as \
                    a comma-separated list): a SAM tag (any pair \
                    of uppercase letters) or mapq, pos5, pos3, \
                    cigar, read_len, matched_bp, algn_ref_span, \
                    algn_read_span, dist_to_5, dist_to_3, seq, \
                    mismatches',
    'max_mismatch':'pairtools dedup:Pairs with both sides mapped within this\
                distance (bp) from each other are considered\
                    duplicates. [dedup option]  [default: 3'
}