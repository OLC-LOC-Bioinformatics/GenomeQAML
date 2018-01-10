# Tests for feature extraction of OLC Quality Assessment Tool
import os

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

import extract_features
from extract_features import *


gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, \
    n75_dict, l50_dict, l75_dict = extract_features.main('tests/test_fastas')


def test_genome_length_normal_case():  # Given a normal(ish) fasta, tests that length is found correctly
    assert genome_length_dict['normal'] == 63


def test_genome_length_with_lowercase():
    assert genome_length_dict['normal_with_lowercase'] == 63


def test_genome_length_blank_contig():
    assert genome_length_dict['blank_contig'] == 55


def test_genome_length_one_contig():
    assert find_genome_length('tests/test_fastas/one_contig.fasta') == 10


def test_find_n50():
    assert find_n50('tests/test_fastas/normal.fasta') == 40


def test_find_n50_with_lowercase():
    assert find_n50('tests/test_fastas/normal_with_lowercase.fasta') == 40


def test_find_n50_blank_contig():
    assert find_n50('tests/test_fastas/blank_contig.fasta') == 40


def test_find_n50_one_contig():
    assert find_n50('tests/test_fastas/one_contig.fasta') == 10


def test_find_n50_lots_of_contigs():
    assert find_n50('tests/test_fastas/several_contigs.fasta') == 7


def test_find_n75():
    assert find_n75('tests/test_fastas/normal.fasta') == 15


def test_find_n75_with_lowercase():
    assert find_n75('tests/test_fastas/normal_with_lowercase.fasta') == 15


def test_find_n75_blank_contig():
    assert find_n75('tests/test_fastas/normal.fasta') == 15


def test_find_n75_one_contig():
    assert find_n75('tests/test_fastas/one_contig.fasta') == 10


def test_find_n75_lots_of_contigs():
    assert find_n75('tests/test_fastas/several_contigs.fasta') == 5


def test_find_num_contigs():
    assert find_num_contigs('tests/test_fastas/normal.fasta') == 3


def test_find_num_contigs_with_lowercase():
    assert find_num_contigs('tests/test_fastas/normal_with_lowercase.fasta') == 3


def test_find_num_contigs_blank_contig():
    assert find_num_contigs('tests/test_fastas/blank_contig.fasta') == 3


def test_find_num_contigs_one_contig():
    assert find_num_contigs('tests/test_fastas/one_contig.fasta') == 1


def test_find_largest_contig():
    assert find_largest_contig('tests/test_fastas/normal.fasta') == 40


def test_find_largest_contig_with_lowercase():
    assert find_largest_contig('tests/test_fastas/normal_with_lowercase.fasta') == 40


def test_find_largest_contig_blank_contig():
    assert find_largest_contig('tests/test_fastas/blank_contig.fasta') == 40


def test_find_largest_contig_one_contig():
    assert find_largest_contig('tests/test_fastas/one_contig.fasta') == 10


def test_find_gc_percent():
    assert find_gc_percent('tests/test_fastas/fifty_gc.fasta') == 50.0


def test_find_gc_percent_twentyfive():
    assert find_gc_percent('tests/test_fastas/twentyfive_gc.fasta') == 25.0


def test_find_gc_percent_seventyfive():
    assert find_gc_percent('tests/test_fastas/twentyfive_gc.fasta') == 75.0


def test_find_l50():
    assert find_l50('tests/test_fastas/normal.fasta') == 1


def test_find_l50_with_lowercase():
    assert find_l50('tests/test_fastas/normal_with_lowercase.fasta') == 1


def test_find_l50_blank_contig():
    assert find_l50('tests/test_fastas/blank_contig.fasta') == 1


def test_find_l50_several_contigs():
    assert find_l50('tests/test_fastas/several_contigs.fasta') == 3


def test_find_l75():
    assert find_l75('tests/test_fastas/normal.fasta') == 2


def test_find_l75_lowercase():
    assert find_l75('tests/test_fastas/normal.fasta') == 2


def test_find_l75_blank_contig():
    assert find_l75('tests/test_fastas/normal.fasta') == 2


def test_find_l75_one_contig():
    assert find_l75('tests/test_fastas/one_contig.fasta') == 1


def test_find_l75_lots_of_contigs():
    assert find_l75('tests/test_fastas/several_contigs.fasta') == 5
