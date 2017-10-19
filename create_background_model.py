import csv
import lib.tf_analysis_tools as tools
import numpy as np
import os

from Bio import SeqIO
from collections import Counter
from multiprocessing import Pool
from pandas import read_csv
from argparse import ArgumentParser
from imp import load_source


CHROMOSOMES = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', \
                'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', \
                'chrX', 'chrY']

parser = ArgumentParser(description='This script analyzes ATAC-Seq and GRO-Seq data and produces various plots for further data analysis.', epilog='IMPORTANT: Please ensure that ALL bed files used with this script are sorted by the same criteria.')
parser.add_argument('-c', '--cell-type', dest='output_prefix', metavar='CELL_TYPE', \
                    help='Cell type (k562, imr90, etc), primarily used for the output file prefix', required=True)
parser.add_argument('-f', '--config', dest='config_filename', metavar='CONFIG_FILE', \
                    help='Configuration file name in ./config, defaults otherwise to "configparams.py"', \
                    default='configparams.py', required=False)
args = parser.parse_args()
cfg = load_source('config.%s' % args.config_filename.split('.')[0], args.config_filename)


# split by chromosomes for multiprocessing
def find_seqs_per_chrom(current_chrom):
    atac_df = read_csv(cfg.atac_overlaps_filename, header=None, sep=",", \
                       usecols=[0, 1, 2, 3], \
                       names=['chrom', 'start', 'end', 'overlaps_tfit'],
                       dtype={'chrom':'str', 'start':'str', 'end':'str', 'overlaps_tfit':'str'})
    atac_iter = atac_df[(atac_df.chrom == current_chrom)].itertuples()
    active_peaks_file = '%s/peaks_%s.bed' % (cfg.fasta_dir, current_chrom)
    output_file = open(active_peaks_file, 'w')
    writer = csv.writer(output_file, delimiter ='\t')
    for atac_peak in atac_iter:
        if atac_peak.overlaps_tfit == '1':
            # create a bedfile of all peaks overlapping TFit
            writer.writerow([atac_peak.chrom, atac_peak.start, atac_peak.end])
    output_file.close()

    # gather the sequence for each of these peaks
    output = os.popen('%s/seqtk subseq %s/%s.fa %s > %s.fa' % (cfg.seqtk_dir, cfg.fasta_dir, current_chrom, active_peaks_file, active_peaks_file)).read()
    fasta_sequences = SeqIO.parse(open('%s.fa' % active_peaks_file),'fasta')
    seqs = []
    for fasta in fasta_sequences:
        if fasta.seq:
            seqs.append(str(fasta.seq))
    print 'done with %s' %  current_chrom
    return seqs


def get_background_per_chrom(current_chrom):
    atac_df = read_csv(cfg.atac_overlaps_filename, header=None, sep=",", \
                       usecols=[0, 1, 2, 3], \
                       names=['chrom', 'start', 'end', 'overlaps_tfit'],
                       dtype={'chrom':'str', 'start':'str', 'end':'str', 'overlaps_tfit':'str'})
    atac_iter = atac_df[(atac_df.chrom == current_chrom)].itertuples()
    background_file = '%s/background_%s.bed' % (cfg.fasta_dir, current_chrom)
    output_file = open(background_file, 'w')
    writer = csv.writer(output_file, delimiter ='\t')
    for atac_peak in atac_iter:
        # create a bedfile of non-peak regions, 100bp-wide
        writer.writerow([atac_peak.chrom, (int(atac_peak.start) - 100), atac_peak.start])
        writer.writerow([atac_peak.chrom, atac_peak.end, (int(atac_peak.end) + 100)])
    output_file.close()

    # gather the sequence for each of these peaks
    output = os.popen('%s/seqtk subseq %s/%s.fa %s > %s.fa' % (cfg.seqtk_dir, cfg.fasta_dir, current_chrom, background_file, background_file)).read()
    fasta_sequences = SeqIO.parse(open('%s.fa' % background_file),'fasta')
    seqs = []
    for fasta in fasta_sequences:
        if fasta.seq:
            seqs.append(str(fasta.seq))
    print 'done with background for %s' %  current_chrom
    return seqs


# A: 0
# C: 1
# T: 2
# G: 3
def get_base_index(nucleotide):
    base = 0        # A
    if nucleotide == 'c':
        base = 1
    elif nucleotide == 't':
        base = 2
    elif nucleotide == 'g':
        base = 3
    #TODO: special symbols? (R, N, etc)
    #      for now, flip a coin on the possible bases the symbol represents
    elif nucleotide == 'r':     # purine
        base = np.random.choice([0,3])
    elif nucleotide == 'y':     # pyrimidine
        base = np.random.choice([1,2])
    elif nucleotide == 'k':     # keto
        base = np.random.choice([2,3])
    elif nucleotide == 'm':     # amino
        base = np.random.choice([0,1])
    elif nucleotide == 's':     # strong
        base = np.random.choice([1,3])
    elif nucleotide == 'w':     # weak
        base = np.random.choice([0,2])
    elif nucleotide == 'b':
        base = np.random.choice([1,2,3])
    elif nucleotide == 'd':
        base = np.random.choice([0,2,3])
    elif nucleotide == 'h':
        base = np.random.choice([0,1,2])
    elif nucleotide == 'v':
        base = np.random.choice([0,1,3])
    elif nucleotide == 'n':     # any
        base = np.random.choice([0,1,2,3])
    return base



if __name__=='__main__':
    all_sequences = find_seqs_per_chrom(CHROMOSOMES[20])
    pool = Pool()
    results = pool.map(find_seqs_per_chrom, CHROMOSOMES)
    pool.close()
    pool.join()
    all_sequences = np.hstack(results)
    print len(all_sequences)
#    all_sequences = [x for x in results if x is not None]

    bg_pool = Pool()
    results = bg_pool.map(get_background_per_chrom, CHROMOSOMES)
    bg_pool.close()
    bg_pool.join()
    all_bg_sequences = np.hstack(results)

    np.save('atac_tfit_overlaps_all_sequences.npy', all_sequences)

# TODO: how do we normalize these, to do a logo? take the longest one and base it on that? (we don't know the polarity


    # transition probability matrix for the mononucleotide distribution
    tpm = np.zeros((4,4), dtype=float)
    bg_tpm = np.zeros((4,4), dtype=float)       # background
    prev_base = np.random.choice([0, 1, 2, 3])
    mu = np.zeros(4, dtype=float)     # initial distribution
    bg_mu = np.zeros(4, dtype=float)     # initial background distribution
    widths = []
    for seq in all_sequences:
        widths.append(len(seq))

        # track the first base of this range to come up with an initial distribution
        first_base = get_base_index(seq[0].lower())
        mu[first_base] += 1

        for nucleotide in seq.lower():
            base = get_base_index(nucleotide)
            tpm[prev_base][base] += 1
            prev_base = base

    for seq in all_bg_sequences:
        # track the first base of this range to come up with an initial distribution
        first_base = get_base_index(seq[0].lower())
        bg_mu[first_base] += 1

        for nucleotide in seq.lower():
            base = get_base_index(nucleotide)
            bg_tpm[prev_base][base] += 1
            prev_base = base
        print tpm

    for i in range(4):
        tpm[i][:] = tpm[i][:] / tpm[i][:].sum()
    mu = mu / mu.sum()

    for i in range(4):
        bg_tpm[i][:] = bg_tpm[i][:] / bg_tpm[i][:].sum()
    bg_mu = bg_mu / bg_mu.sum()

    np.save('atac_tfit_overlaps_tpm.npy', tpm)
    np.save('atac_tfit_overlaps_initial_distribution.npy', mu)

    np.save('atac_background_tpm.npy', bg_tpm)
    np.save('atac_background_initial_distribution.npy', bg_mu)

    print 'Transition probability matrix for the mononucleotide model:'
    print tpm

    mean_width = np.mean(widths)
    std_width = np.std(widths)
    total_peak_count = os.popen('wc -l %s' % cfg.atac_overlaps_filename).read().split()[0]
    np.save('atac_peak_stats.npy', [total_peak_count, mean_width, std_width])
    print 'ATAC peaks width ~ Normal(%.3f, %.3f)' % (mean_width, std_width)



# figure out dinucleotide distribution



