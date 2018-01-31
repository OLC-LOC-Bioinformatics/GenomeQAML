# Tests for feature extraction of OLC Quality Assessment Tool
import os

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

from genomeqaml import extract_features


gc_dict, contig_dist_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, n75_dict, \
    n90_dict, l50_dict, l75_dict, l90_dict, orf_dist_dict = extract_features.main('tests/test_fastas',
                                                                                  refseq_database='refseq.msh',
                                                                                  report=False)


def test_genome_length_normal_case():  # Given a normal(ish), tests that length is found correctly
    assert genome_length_dict['normal'] == 63


def test_genome_length_with_lowercase():
    assert genome_length_dict['normal_with_lowercase'] == 63


def test_genome_length_blank_contig():
    assert genome_length_dict['blank_contig'] == 55


def test_genome_length_one_contig():
    assert genome_length_dict['one_contig'] == 10


def test_find_n50():
    assert n50_dict['normal'] == 40


def test_find_n50_with_lowercase():
    assert n50_dict['normal_with_lowercase'] == 40


def test_find_n50_blank_contig():
    assert n50_dict['blank_contig'] == 40


def test_find_n50_one_contig():
    assert n50_dict['one_contig'] == 10


def test_find_n50_lots_of_contigs():
    assert n50_dict['several_contigs'] == 7


def test_find_n75():
    assert n75_dict['normal'] == 15


def test_find_n75_with_lowercase():
    assert n75_dict['normal_with_lowercase'] == 15


def test_find_n75_blank_contig():
    assert n75_dict['normal'] == 15


def test_find_n75_one_contig():
    assert n75_dict['one_contig'] == 10


def test_find_n75_lots_of_contigs():
    assert n75_dict['several_contigs'] == 5


def test_num_contigs_dict():
    assert num_contigs_dict['normal'] == 3


def test_num_contigs_dict_with_lowercase():
    assert num_contigs_dict['normal_with_lowercase'] == 3


def test_num_contigs_dict_blank_contig():
    assert num_contigs_dict['blank_contig'] == 3


def test_num_contigs_dict_one_contig():
    assert num_contigs_dict['one_contig'] == 1


def test_longest_contig_dict():
    assert longest_contig_dict['normal'] == 40


def test_longest_contig_dict_with_lowercase():
    assert longest_contig_dict['normal_with_lowercase'] == 40


def test_longest_contig_dict_blank_contig():
    assert longest_contig_dict['blank_contig'] == 40


def test_longest_contig_dict_one_contig():
    assert longest_contig_dict['one_contig'] == 10


def test_gc_dict():
    assert gc_dict['fifty_gc'] == 50.0


def test_gc_dict_twentyfive():
    assert gc_dict['twentyfive_gc'] == 25.0


def test_gc_dict_seventyfive():
    assert gc_dict['seventyfive_gc'] == 75.0


def test_l50_dict():
    assert l50_dict['normal'] == 1


def test_l50_dict_with_lowercase():
    assert l50_dict['normal_with_lowercase'] == 1


def test_l50_dict_blank_contig():
    assert l50_dict['blank_contig'] == 1


def test_l50_dict_several_contigs():
    assert l50_dict['several_contigs'] == 3


def test_l75_dict():
    assert l75_dict['normal'] == 2


def test_l75_dict_lowercase():
    assert l75_dict['normal_with_lowercase'] == 2


def test_l75_dict_blank_contig():
    assert l75_dict['normal'] == 2


def test_l75_dict_one_contig():
    assert l75_dict['one_contig'] == 1


def test_l75_dict_lots_of_contigs():
    assert l75_dict['several_contigs'] == 5


def test_l90_dict_normal():
    assert l90_dict['normal'] == 3


def test_l90_dict_normal_with_lowercase():
    assert l90_dict['normal_with_lowercase'] == 3


def test_l90_dict_blank_contig():
    assert l90_dict['blank_contig'] == 2


def test_n90_dict_normal():
    assert n90_dict['normal'] == 8
