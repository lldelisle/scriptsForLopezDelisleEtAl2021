import numpy as np
from scipy.stats import lognorm
import argparse
import sys


def generate_data(ncells, scale,
                  mean,
                  seed, fo):
    data = lognorm.rvs(s=scale, scale=mean, size=ncells, random_state=seed)
    np.savetxt(fo, data, fmt="%d", header='nCount_RNA', comments='')


def parse_arguments(args=None):
    argp = argparse.ArgumentParser(
        description=("Generate counts following a log-normal distribution."))
    argp.add_argument('--ncells', default=None, required=True,
                      type=int,
                      help="Number of cells.")
    argp.add_argument('--scale', default=0.3, type=float,
                      help="Scale of the log-normal distribution.")
    argp.add_argument('--mean', default=16000, type=float,
                      help="Average of the log-normal distribution.")
    argp.add_argument('--seed', default=1, type=int,
                      help="Set seed for reproductible results.")
    argp.add_argument('--output', default=sys.stdout,
                      type=argparse.FileType('w'),
                      help="Output table.")
    return(argp)


def main(args=None):
    args = parse_arguments().parse_args(args)
    generate_data(args.ncells, args.scale,
                  args.mean,
                  args.seed, args.output)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
