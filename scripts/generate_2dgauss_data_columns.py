import numpy as np
from scipy.stats import poisson, uniform
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


def trunc_norm_2d(mu, sigma, corr, size, seed):
  try:
    rng = np.random.default_rng(seed)
  except AttributeError:
    # For older numpy versions:
    np.random.seed(seed)
    rng = np.random
  cov = np.array([[sigma[0] * sigma[0], sigma[0] * sigma[1] * corr],
                  [sigma[0] * sigma[1] * corr, sigma[1] * sigma[1]]])
  values = rng.multivariate_normal(mu, cov, size)
  mask_0 = [v[0] < 0 or v[1] < 0 for v in values]
  # Because we want only positive expression:
  while sum(mask_0) > 0:
    values[mask_0] = rng.multivariate_normal(mu, cov, sum(mask_0))
    mask_0 = [v[0] < 0 or v[1] < 0 for v in values]
  return(values)


def generate_values_from_name(colname, n, seed):
  cn_split = colname.split('_')
  amp = np.array([float(p) for p in cn_split[::6]])
  # Amplitude
  if np.any(amp < 0) or np.sum(amp) > 1:
    raise Exception("Negative amplitude")
  mu = np.array([[float(mux), float(muy)] for mux, muy in zip(cn_split[1::6], cn_split[2::6])])
  sigma = np.array([[float(sigx), float(sigy)] for sigx, sigy in zip(cn_split[3::6], cn_split[4::6])])
  corr = np.array([float(p) for p in cn_split[5::6]]) 
  masks = get_masks(amp, n, seed)
  exp_values = np.zeros((n, 2))
  for i, (cur_mu, cur_sigma, cur_corr) in enumerate(zip(mu, sigma, corr)):
    exp_values[masks[i]] = trunc_norm_2d(mu=cur_mu, sigma=cur_sigma,
                                         corr=cur_corr,
                                         size=sum(masks[i]),
                                         seed=seed)
  return(exp_values)


def generate_2dgauss_data(input_file,
                          colnames,
                          starting_seed, fo):
    data = pd.read_csv(input_file, sep="\t")
    N = data['nCount_RNA']
    for id_col, colname in enumerate(colnames):
      exp_values = generate_values_from_name(colname, N.size, starting_seed + id_col)
      exp_values_x, exp_values_y  = np.transpose(exp_values)
      # We now apply the poisson distribution
      ks_x = poisson.rvs(mu=N * 1e-4 * (np.exp(exp_values_x) - 1),
                         random_state=starting_seed + id_col)
      ks_y = poisson.rvs(mu=N * 1e-4 * (np.exp(exp_values_y) - 1),
                         random_state=starting_seed + id_col)
      data[f'{colname}_x'] = ks_x
      data[f'{colname}_y'] = ks_y
      # We also store the expression before poisson
      data[f'{colname}_x_expression'] = exp_values_x
      data[f'{colname}_y_expression'] = exp_values_y
    # Export
    data.to_csv(fo, index = False, header=True, sep='\t')


def parse_arguments(args=None):
    argp = argparse.ArgumentParser(
        description=("Generate values x and y following known distributions composed of"
                     " a prop of zeros + 2d normal distribution(s)."))
    argp.add_argument('--input', default=None, required=True,
                      help="Path for a table with at least one column called 'nCount_RNA'.")
    argp.add_argument('--colnames', default=None, required=True, nargs='+',
                      help="colnames to simulate (format is "
                      "{amp1}_{mux1}_{muy1}_{sigx1}_{sigy1}[_{amp2}_{mux2}_{muy2}_{sigx2}_{sigy2}]...)."
                      " For gauss, mu and sig are the mean and the sigma of the "
                      "non-trucated gaussian. Amp is the proportion of cells in this distribution."
                      " The proportion of 0 is 1 - sum(amp1, amp2...)")
    argp.add_argument('--startingSeed', default=1, type=int,
                      help="Set seed for reproductible results "
                      "(for each column the seed used for masks is startingSeed + 5 * col_id).")
    argp.add_argument('--output', default=sys.stdout,
                      type=argparse.FileType('w'),
                      help="Output table.")
    return(argp)


def main(args=None):
    args = parse_arguments().parse_args(args)
    generate_2dgauss_data(args.input, args.colnames,
                          args.startingSeed, args.output)

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
