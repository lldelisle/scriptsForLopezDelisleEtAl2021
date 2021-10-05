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


def generate_values_from_name(colname, n, seed, xscale):
  distrib = colname.split('_')[::4]
  amp = np.array([float(p) for p in colname.split('_')[1::4]])
  # Amplitude
  if np.any(amp < 0) or np.sum(amp) > 1:
    raise Exception("Negative amplitude")
  loc = np.array([float(p) for p in colname.split('_')[2::4]])
  scale = np.array([float(p) for p in colname.split('_')[3::4]])
  masks = get_masks(amp, n, seed)
  exp_values = np.zeros(n)
  if xscale == 'log':
    exp_values[:] = -np.Inf
  for i, (cur_distrib, cur_loc, cur_scale) in enumerate(zip(distrib, loc, scale)):
    if cur_distrib == "uniform":
      exp_values[masks[i]] = uniform.rvs(loc=cur_loc, scale=cur_scale,
                                         size=sum(masks[i]),
                                         random_state=seed)
    elif cur_distrib == "gauss":
          if xscale == 'Seurat':
            exp_values[masks[i]] = truncnorm.rvs(- cur_loc / cur_scale, np.inf,
                                                 loc=cur_loc, scale=cur_scale,
                                                 size=sum(masks[i]),
                                                 random_state=seed)
          elif xscale == 'log':
            exp_values[masks[i]] = norm.rvs(loc=cur_loc, scale=cur_scale,
                                            size=sum(masks[i]),
                                            random_state=seed)
          else:
            raise Exception("Only Seurat and log are implemented.")
    else:
      raise Exception("The distribution asked is not implemented.")
  return(exp_values)


def generate_paired_data(input_file,
                         colname_x1,
                         colname_x2,
                         colname_y1,
                         colname_y2,
                         props4groups,
                         starting_seed, fo,
                         xscale, outputSanity):
    data = pd.read_csv(input_file, sep="\t")
    N = data['nCount_RNA']
    # There are 4 groups:
    # g0: x1  y1
    # g1: x1  y2
    # g2: x2  y1
    # g3: x1  y2
  
    # We generate data whose final distribution would be:
    # x = colname_x1 + colname_x2
    # y = colname_y1 + colname_y2

    for id_col, prop4group in enumerate(props4groups):
      try:
        amps = [float(amp) for amp in prop4group.split('_')]
        assert sum(amps) == 1, "The sum of proportions should be 1."
      except Exception as e:
        sys.stderr.write(f"The prop4group {prop4group} could not be computed: {e}")
        continue
      masks = get_masks(amps, N.size, starting_seed + 5 * id_col)
      exp_values_x = np.zeros(N.size)
      exp_values_y = np.zeros(N.size)
      # We don't need to set to -np.Inf for log
      # Because we are sure that the masks cover all.
      # Generate x for g0, g1
      exp_values_x[masks[0] + masks[1]] = \
        generate_values_from_name(colname_x1, sum(masks[0] + masks[1]), starting_seed + 5 * id_col + 1, xscale)
      # Generate x for g2, g3
      exp_values_x[masks[2] + masks[3]] = \
        generate_values_from_name(colname_x2, sum(masks[2] + masks[3]), starting_seed + 5 * id_col + 2, xscale)
      # Generate y for g0, g2
      exp_values_y[masks[0] + masks[2]] = \
        generate_values_from_name(colname_y1, sum(masks[0] + masks[2]), starting_seed + 5 * id_col + 3, xscale)
      # Generate y for g1, g3
      exp_values_y[masks[1] + masks[3]] = \
        generate_values_from_name(colname_y2, sum(masks[1] + masks[3]), starting_seed + 5 * id_col + 4, xscale)
      # We now apply the poisson distribution
      if xscale == 'Seurat':
        ks_x = poisson.rvs(mu=N * 1e-4 * (np.exp(exp_values_x) - 1),
                           random_state=starting_seed + id_col)
        ks_y = poisson.rvs(mu=N * 1e-4 * (np.exp(exp_values_y) - 1),
                          random_state=starting_seed + id_col)
      elif xscale == 'log':
        ks_x = poisson.rvs(mu=N * np.exp(exp_values_x),
                           random_state=starting_seed + id_col)
        ks_y = poisson.rvs(mu=N * np.exp(exp_values_y),
                           random_state=starting_seed + id_col)
      else:
        raise Exception("Only Seurat and log are implemented.")
      data[f'{prop4group}_x'] = ks_x
      data[f'{prop4group}_y'] = ks_y
      # We also store the expression before poisson
      data[f'{prop4group}_x_expression'] = exp_values_x
      data[f'{prop4group}_y_expression'] = exp_values_y
      # We also store the group id:
      group = np.zeros(N.size)
      group[masks[1]] = 1
      group[masks[2]] = 2
      group[masks[3]] = 3
      data[f'{prop4group}_group'] = group
    # Export
    data.to_csv(fo, index = False, header=True, sep='\t')

    # Export for Sanity
    if outputSanity is not None:
      count_columns = [a + b for a, b in zip(np.repeat(props4groups, 2), ['_x', '_y'] * len(props4groups))]
      data['complement'] = data.loc[:, 'nCount_RNA'] - np.sum(data.loc[:, count_columns], axis=1)
      data.loc[:, count_columns + ['complement']].T.to_csv(outputSanity, index=True, header=True, sep='\t', index_label = 'GeneID')


def parse_arguments(args=None):
    argp = argparse.ArgumentParser(
        description=("Generate values x and y following known distributions composed of"
                     " a prop of zeros + uniform or normal distribution(s)."
                     " x is the sum of 2 distributions, y also. Cells are split in"
                     " 4 groups corresponding to the 4 combinations of distibutions"
                     " for x and y."))
    argp.add_argument('--input', default=None, required=True,
                      help="Path for a table with at least one column called 'nCount_RNA'.")
    argp.add_argument('--colnamex1', default=None, required=True,
                      help="colname to simulate for x group 0 and 1(format is "
                      "[uniform|gauss]_{amp1}_{loc1}_{scale1}[_{amp2}_{loc2}_{scale2}]...)."
                      " For gauss, loc and scale are the mean and the sigma of the "
                      "non-trucated gaussian. For uniform, loc is the min and scale"
                      " the length. Amp is the proportion of cells in this distribution."
                      " The proportion of 0 is 1 - sum(amp1, amp2...)")
    argp.add_argument('--colnamex2', default=None, required=True,
                      help="colname to simulate for x group 2 and 3(format is "
                      "[uniform|gauss]_{amp1}_{loc1}_{scale1}[_{amp2}_{loc2}_{scale2}]...)."
                      " For gauss, loc and scale are the mean and the sigma of the "
                      "non-trucated gaussian. For uniform, loc is the min and scale"
                      " the length. Amp is the proportion of cells in this distribution."
                      " The proportion of 0 is 1 - sum(amp1, amp2...)")
    argp.add_argument('--colnamey1', default=None, required=True,
                      help="colname to simulate for x group 0 and 2(format is "
                      "[uniform|gauss]_{amp1}_{loc1}_{scale1}[_{amp2}_{loc2}_{scale2}]...)."
                      " For gauss, loc and scale are the mean and the sigma of the "
                      "non-trucated gaussian. For uniform, loc is the min and scale"
                      " the length. Amp is the proportion of cells in this distribution."
                      " The proportion of 0 is 1 - sum(amp1, amp2...)")
    argp.add_argument('--colnamey2', default=None, required=True,
                      help="colname to simulate for x group 1 and 3(format is "
                      "[uniform|gauss]_{amp1}_{loc1}_{scale1}[_{amp2}_{loc2}_{scale2}]...)."
                      " For gauss, loc and scale are the mean and the sigma of the "
                      "non-trucated gaussian. For uniform, loc is the min and scale"
                      " the length. Amp is the proportion of cells in this distribution."
                      " The proportion of 0 is 1 - sum(amp1, amp2...)")
    argp.add_argument('--props4groups', default=['0.25_0.25_0.25_0.25'], required=True, nargs='+',
                      help="proportions of 4 groups to simulate (format is "
                      "group0_group1_group2_group4. The sum should be 1.")
    argp.add_argument('--xscale', default='Seurat', choices=['Seurat', 'log'],
                      help="Scale used to convert expression to lambda in Poisson sampling.")
    argp.add_argument('--startingSeed', default=1, type=int,
                      help="Set seed for reproductible results "
                      "(for each column the seed used for masks is startingSeed + 5 * col_id).")
    argp.add_argument('--output', default=sys.stdout,
                      type=argparse.FileType('w'),
                      help="Output table.")
    argp.add_argument('--outputSanity', default=None,
                      help="Output table for Sanity.")
    return(argp)


def main(args=None):
    args = parse_arguments().parse_args(args)
    generate_paired_data(args.input, args.colnamex1, args.colnamex2, 
                         args.colnamey1, args.colnamey2,
                         args.props4groups, args.startingSeed, args.output,
                         args.xscale, args.outputSanity)

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
