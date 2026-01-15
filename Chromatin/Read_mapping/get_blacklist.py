import sys
import math
import argparse
import pyBigWig


def get_args(argv):
    """
    Parses arguments from command line
    """
    args = argparse.ArgumentParser()
    args.add_argument('-i', required=True,
                      metavar='<bigWig>',
                      help='input bigWig file')
    args.add_argument('-t', default=10,
                      type = float,
                      metavar='<num>',
                      help='cutoff as a fold SD over the mean (default=10)')
    args.add_argument('-o', required=True,
                      metavar='<bed>',
                      help='output bed file')
    return args.parse_args()


def high_dens(bw, cutoff, outfile):
    """
    Writes out regions with high density of reads
    """
    bw_mean = bw.header()['sumData'] / bw.header()['nBasesCovered']
    bw_std = math.sqrt(bw.header()['sumSquared']/bw.header()['nBasesCovered'] - bw_mean**2)
    print(bw_mean)
    print(bw_std)
    cut_value = bw_mean + bw_std * cutoff
    out = open(outfile,'w')
    for chrom in bw.chroms().keys():
        intervals = bw.intervals(chrom)
        start = -1
        end = -1

        for i in intervals:
            if i[2] > cut_value:
                if i[0] == end:
                    end = i[1]
                else:
                    if start > -1:
                        out.write(chrom + '\t' + str(start) + '\t' + str(end) + '\n')
                    start = i[0]
                    end = i[1]       
    out.close()


if __name__ == '__main__':
    args = get_args(sys.argv[1:])
    bw = pyBigWig.open(args.i)
    high_dens(bw, args.t, args.o)
