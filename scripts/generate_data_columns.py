import numpy as np
from scipy.stats import poisson, uniform, truncnorm, norm
import argparse
import sys
import pandas as pd
import random


def get_masks(amp, n, seed):
  random = uniform.rvs(size=n, random_state=seed)
  masks = []
  cur_threshold = 0
  for cur_amp in amp:
    masks.append((random >= cur_threshold) & (random < cur_threshold + cur_amp))
    cur_threshold += cur_amp
  masks.append(random >= cur_threshold)
  return(masks)


def generate_data(input_file, colnames, starting_seed,
                  fo, xscale, outputSanity):
    data = pd.read_csv(input_file, sep="\t")
    N = data['nCount_RNA']
    for id_col, colname in enumerate(colnames):
      distrib = colname.split('_')[::4]
      amp = np.array([float(p) for p in colname.split('_')[1::4]])
      # Amplitude
      if np.any(amp < 0) or np.sum(amp) > 1:
        raise Exception("Negative amplitude")
      if xscale == 'log' and np.sum(amp) != 1:
        sys.stderr.write("Warning: When using log scale 0 become -Inf.\n")
      loc = np.array([float(p) for p in colname.split('_')[2::4]])
      scale = np.array([float(p) for p in colname.split('_')[3::4]])
      masks = get_masks(amp, N.size, starting_seed + id_col)
      exp_values = np.zeros(N.size)
      if xscale == 'log':
        exp_values[:] = -np.Inf
      for i, (cur_distrib, cur_loc, cur_scale) in enumerate(zip(distrib, loc, scale)):
        if cur_distrib == "uniform":
          exp_values[masks[i]] = uniform.rvs(loc=cur_loc, scale=cur_scale,
                                             size=sum(masks[i]),
                                             random_state=starting_seed + id_col)
        elif cur_distrib == "gauss":
          if xscale == 'Seurat':
            exp_values[masks[i]] = truncnorm.rvs(- cur_loc / cur_scale, np.inf,
                                                 loc=cur_loc, scale=cur_scale,
                                                 size=sum(masks[i]),
                                                 random_state=starting_seed + id_col)
          elif xscale == 'log':
            exp_values[masks[i]] = norm.rvs(loc=cur_loc, scale=cur_scale,
                                            size=sum(masks[i]),
                                            random_state=starting_seed + id_col)
          else:
            raise Exception("Only Seurat and log are implemented.")
        else:
          raise Exception("The distribution asked is not implemented.")
      # We now apply the poisson distribution
      if xscale == 'Seurat':
        ks = poisson.rvs(mu=N * 1e-4 * (np.exp(exp_values) - 1),
                         random_state=starting_seed + id_col)
      elif xscale == 'log':
        ks = poisson.rvs(mu=N * np.exp(exp_values),
                         random_state=starting_seed + id_col)
      else:
        raise Exception("Only Seurat and log are implemented.")
      data[colname] = ks
      # We also store the expression before poisson
      data[f'{colname}_expression'] = exp_values

    # Export
    data.to_csv(fo, index = False, header=True, sep='\t')

    # Export for Sanity
    if outputSanity is not None:
      data['complement'] = data.loc[:, 'nCount_RNA'] - np.sum(data.loc[:, colnames], axis=1)
      data.loc[:, colnames + ['complement']].T.to_csv(outputSanity, index=True, header=True, sep='\t', index_label = 'GeneID')

def parse_arguments(args=None):
    argp = argparse.ArgumentParser(
        description=("Generate values following a known distribution composed of"
                     " a prop of zeros + uniform or normal distribution(s)."))
    argp.add_argument('--input', default=None, required=True,
                      help="Path for a table with at least one column called 'nCount_RNA'.")
    argp.add_argument('--colnames', default=None, required=True, nargs='+',
                      help="colnames to simulate (format is "
                      "[uniform|gauss]_{amp1}_{loc1}_{scale1}[_{amp2}_{loc2}_{scale2}]...)."
                      " For gauss, loc and scale are the mean and the sigma of the "
                      "non-trucated gaussian. For uniform, loc is the min and scale"
                      " the length. Amp is the proportion of cells in this distribution."
                      " The proportion of 0 is 1 - sum(amp1, amp2...)")
    argp.add_argument('--xscale', default='Seurat', choices=['Seurat', 'log'],
                      help="Scale used to convert expression to lambda in Poisson sampling.")
    argp.add_argument('--startingSeed', default=1, type=int,
                      help="Set seed for reproductible results "
                      "(for each column the seed used is startingSeed + col_id).")
    argp.add_argument('--output', default=sys.stdout,
                      type=argparse.FileType('w'),
                      help="Output table.")
    argp.add_argument('--outputSanity', default=None,
                      help="Output table for Sanity.")
    return(argp)


def main(args=None):
    args = parse_arguments().parse_args(args)
    generate_data(args.input, args.colnames,
                  args.startingSeed, args.output,
                  args.xscale, args.outputSanity)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
