PAIRSAM_FORMAT_VERSION = "1.0.0"

PAIRSAM_SEP = "\t"
PAIRSAM_SEP_ESCAPE = r"\t"
SAM_SEP = "\031"
SAM_SEP_ESCAPE = r"\031"
INTER_SAM_SEP = "\031NEXT_SAM\031"

COL_READID = 0
COL_C1 = 1
COL_P1 = 2
COL_C2 = 3
COL_P2 = 4
COL_S1 = 5
COL_S2 = 6
COL_PTYPE = 7
COL_SAM1 = 8
COL_SAM2 = 9

COLUMNS = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
    "sam1",
    "sam2",
    "walk_pair_index",
    "walk_pair_type",
]

# Required columns for formats:
COLUMNS_PAIRSAM = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
    "sam1",
    "sam2",
]

COLUMNS_PAIRS = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
]

UNMAPPED_CHROM = "!"
UNMAPPED_POS = 0
UNMAPPED_STRAND = "-"

UNANNOTATED_RFRAG = -1

EXTRA_COLUMNS = [
    "mapq",
    "pos5",
    "pos3",
    "cigar",
    "read_len",
    "matched_bp",
    "algn_ref_span",
    "algn_read_span",
    "dist_to_5",
    "dist_to_3",
    "seq",
    "mismatches",  # Format: "{ref_letter}:{mut_letter}:{phred}:{ref_position}:{read_position}"
]
