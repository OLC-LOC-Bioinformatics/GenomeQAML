#!/usr/bin/env python 3
from Bio.SeqUtils import GC
from Bio import SeqIO
from glob import glob
import click
import os
__author__ = 'adamkoziol'


def main(sequencepath):
    """
    Run the appropriate functions in order
    :param sequencepath: Path of folder containing FASTA genomes
    :return: gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, n75_dict, l50_dict, l75_dict
    """
    files = find_files(sequencepath)
    file_dict = filer(files)
    file_records = fasta_records(file_dict)
    contig_len_dict, gc_dict = fasta_stats(file_dict, file_records)
    longest_contig_dict = find_largest_contig(contig_len_dict)
    genome_length_dict = find_genome_length(contig_len_dict)
    num_contigs_dict = find_num_contigs(contig_len_dict)
    n50_dict = find_n50(contig_len_dict, genome_length_dict)
    n75_dict = find_n75(contig_len_dict, genome_length_dict)
    l50_dict = find_l50(contig_len_dict, genome_length_dict)
    l75_dict = find_l75(contig_len_dict, genome_length_dict)
    report(gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, n75_dict, l50_dict, l75_dict,
           sequencepath)
    return gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, n75_dict, l50_dict, l75_dict


def find_files(sequencepath):
    """
    Use glob to find all FASTA files in the provided sequence path. NOTE: FASTA files must have an extension such as
    .fasta, .fa, or .fas. Extensions of .fsa, .tfa, ect. are not currently supported
    :param sequencepath:
    :return: list of FASTA files
    """
    # Create a sorted list of all the FASTA files in the sequence path
    files = sorted(glob(os.path.join(sequencepath, '*.fa*')))
    return files


def filer(filelist):
    """
    Helper script that creates a dictionary of the stain name: /sequencepath/strain_name.extension)
    :param filelist: List of files to parse
    :return filedict: dictionary of stain name: /sequencepath/strain_name.extension
    """
    # Initialise the dictionary
    filedict = dict()
    for seqfile in filelist:
        # Split off the file extension and remove the path from the name
        strainname = os.path.splitext(os.path.basename(seqfile))[0]
        # Populate the dictionary
        filedict[strainname] = seqfile
    return filedict


def fasta_records(files):
    """
    Use SeqIO to create dictionaries of all records for each FASTA file
    :param files: dictionary of stain name: /sequencepath/strain_name.extension
    :return: file_records: dictionary of all contig records for all strains
    """
    # Initialise the dictionary
    file_records = dict()
    for file_name, fasta in files.items():
        # Create a dictionary of records for each file
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        # Set the records dictionary as the value for file_records
        file_records[file_name] = record_dict
    return file_records


def fasta_stats(files, records):
    """
    Parse the lengths of all contigs for each sample, as well as the total GC%
    :param files: dictionary of stain name: /sequencepath/strain_name.extension
    :param records: Dictionary of strain name: SeqIO records
    :return: contig_len_dict, gc_dict: dictionaries of list of all contig length, and total GC% for all strains
    """
    # Initialise dictionaries
    contig_len_dict = dict()
    gc_dict = dict()
    for file_name in files:
        # Initialise variables to store appropriate values parsed from contig records
        contig_lengths = list()
        fasta_sequence = str()
        for contig, record in records[file_name].items():
            # Append the length of the contig to the list
            contig_lengths.append(len(record.seq))
            # Add the contig sequence to the string
            fasta_sequence += record.seq
        # Set the reverse sorted (e.g. largest to smallest) list of contig sizes as the value
        contig_len_dict[file_name] = sorted(contig_lengths, reverse=True)
        # Calculate the GC% of the total genome sequence using GC - format to have two decimal places
        gc_dict[file_name] = float('{:0.2f}'.format(GC(fasta_sequence)))
    return contig_len_dict, gc_dict


def find_largest_contig(contig_lengths_dict):
    """
    Determine the largest contig for each strain
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :return: longest_contig_dict: dictionary of strain name: longest contig
    """
    # Initialise the dictionary
    longest_contig_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        # As the list is sorted in descending order, the largest contig is the first entry in the list
        longest_contig_dict[file_name] = contig_lengths[0]
    return longest_contig_dict


def find_genome_length(contig_lengths_dict):
    """
    Determine the total length of all the contigs for each strain
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :return: genome_length_dict: dictionary of strain name: total genome length
    """
    # Initialise the dictionary
    genome_length_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        # Use the sum() method to add all the contig lengths in the list
        genome_length_dict[file_name] = sum(contig_lengths)
    return genome_length_dict


def find_num_contigs(contig_lengths_dict):
    """
    Count the total number of contigs for each strain
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :return: num_contigs_dict: dictionary of strain name: total number of contigs
    """
    # Initialise the dictionary
    num_contigs_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        # Use the len() method to count the number of entries in the list
        num_contigs_dict[file_name] = len(contig_lengths)
    return num_contigs_dict


def find_n50(contig_lengths_dict, genome_length_dict):
    """
    Calculate the N50 for each strain. N50 is defined as the largest contig such that at least half of the total
    genome size is contained in contigs equal to or larger than this contig
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :param genome_length_dict: dictionary of strain name: total genome length
    :return: n50_dict: dictionary of strain name: N50
    """
    # Initialise the dictionary
    n50_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        # Initialise a variable to store a running total of contig lengths
        currentlength = 0
        for contig_length in contig_lengths:
            # Increment the current length with the length of the current contig
            currentlength += contig_length
            # If the current length is now greater than the total genome / 2, the current contig length is the N50
            if currentlength >= genome_length_dict[file_name] * 0.5:
                # Populate the dictionary, and break the loop
                n50_dict[file_name] = contig_length
                break
    return n50_dict


def find_n75(contig_lengths_dict, genome_length_dict):
    """
    Calculate the N75 for each strain. N75 is defined as the largest contig such that at least 3/4 of the total
    genome size is contained in contigs equal to or larger than this contig
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :param genome_length_dict: dictionary of strain name: total genome length
    :return: n75_dict: dictionary of strain name: N75
    """
    # Initialise the dictionary
    n75_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        currentlength = 0
        for contig_length in contig_lengths:
            currentlength += contig_length
            # If the current length is now greater than the 3/4 of the total genome length, the current contig length
            # is the N75
            if currentlength >= genome_length_dict[file_name] * 0.75:
                n75_dict[file_name] = contig_length
                break
    return n75_dict


def find_l50(contig_lengths_dict, genome_length_dict):
    """
    Calculate the L50 for each strain. L50 is defined as the number of contigs required to achieve the N50
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :param genome_length_dict: dictionary of strain name: total genome length
    :return: l50_dict: dictionary of strain name: L50
    """
    # Initialise the dictionary
    l50_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        currentlength = 0
        # Initialise a variable to count how many contigs have been added to the currentlength variable
        currentcontig = 0
        for contig_length in contig_lengths:
            currentlength += contig_length
            # Increment :currentcontig each time a contig is added to the current length
            currentcontig += 1
            # Same logic as with the N50, but the contig number is added instead of the length of the contig
            if currentlength >= genome_length_dict[file_name] * 0.5:
                l50_dict[file_name] = currentcontig
                break
    return l50_dict


def find_l75(contig_lengths_dict, genome_length_dict):
    """
    Calculate the L50 for each strain. L75 is defined as the number of contigs required to achieve the N75
    :param contig_lengths_dict: dictionary of strain name: reverse-sorted list of all contig lengths
    :param genome_length_dict: dictionary of strain name: total genome length
    :return: l50_dict: dictionary of strain name: L75
    """
    # Initialise the dictionary
    l75_dict = dict()
    for file_name, contig_lengths in contig_lengths_dict.items():
        currentlength = 0
        currentcontig = 0
        for contig_length in contig_lengths:
            currentlength += contig_length
            currentcontig += 1
            # Same logic as with the L75, but the contig number is added instead of the length of the contig
            if currentlength >= genome_length_dict[file_name] * 0.75:
                l75_dict[file_name] = currentcontig
                break
    return l75_dict


def report(gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, n75_dict, l50_dict, l75_dict,
           sequencepath):
    """
    Create a report of all the extracted features
    :param gc_dict: dictionary of strain name: GC%
    :param longest_contig_dict: dictionary of strain name: longest contig
    :param genome_length_dict: dictionary of strain name: total genome length
    :param num_contigs_dict: dictionary of strain name: total number of contigs
    :param n50_dict: dictionary of strain name: N50
    :param n75_dict: dictionary of strain name: N75
    :param l50_dict: dictionary of strain name: L50
    :param l75_dict: dictionary of strain name: L75
    :param sequencepath: path of folder containing FASTA genomes
    """
    # Initialise string with header information
    data = 'SampleName,TotalLength,NumContigs,LongestContig,N50,N75,L50,L75,GC%\n'
    # Create and open the report for writign
    with open(os.path.join(sequencepath, 'extracted_features.csv'), 'w') as feature_report:
        for file_name in longest_contig_dict:
            # Populate the data string with the appropriate values
            data += '{name},{totlen},{numcontigs},{longestcontig},{n50},{n75},{l50},{l75},{gc}\n'\
                .format(name=file_name,
                        totlen=genome_length_dict[file_name],
                        numcontigs=num_contigs_dict[file_name],
                        longestcontig=longest_contig_dict[file_name],
                        n50=n50_dict[file_name],
                        n75=n75_dict[file_name],
                        l50=l50_dict[file_name],
                        l75=l75_dict[file_name],
                        gc=gc_dict[file_name])
        # Write the string to file
        feature_report.write(data)


# Initialise the click decorator
@click.command()
@click.option('-s', '--sequencepath',
              type=click.Path(exists=True),
              required=True,
              help='Path of folder containing multi-FASTA files')
def cli(sequencepath):
    """
    Pass command line arguments to, and run the feature extraction functions
    """
    main(sequencepath)

