# HIPNEX

Matlab implementation of the algorithm proposed in

Alves, M. M.,  J. M. Pereira, and B. F. Svaiter,
[**A search-free $O (1/k^{3/2})$ homotopy inexact
proximal-Newton extragradient algorithm for
monotone variational inequalities**](https://arxiv.org/abs/2308.05887), 
arXiv:2308.05887 (2023).

## Usage

This repository is self contained. After downloading, you may use the function `hipnex` right away to solve unconstrained variational inequality problems. Alternatively, we provide the function `plain_npe` for an implementation of the Newton Proximal-Extragradient method.
See the documentation of `hipnex` and `plain_npe` regarding usage.

### Reproducing results in the paper

To reproduce the experiments in the paper, run `run_experiments` to run the algorithms and then run `plot_experiments` and `table_experiments` to obtain the corresponding plot and table.

### Other examples

Besides the experiments included in the paper, we also include in this repository a quadratic min-max experiment, including its setup, so that users can more easily understand how to set-up an optimization problem to use HIPNEX or NPE.

## External software

Although the repository is self-contained, we include software from other authors, which we acknowledge here. Namely, we include: 
  - `ORN_ls_simple`: code for Second Order Generalized Optimistic method, proposed in:  
  R. Jiang, A. Mokhtari, and J. M. Pereira, and B. F. Svaiter,
  [**Generalized Optimistic methods for convex-concave saddle point Problems**](https://arxiv.org/abs/2202.09674), arXiv:2202.09674 (2022).
  We thank [Ruichen Jiang](https://github.com/Raymond30/) for providing the code after personal request, and allowing it to be shared in this repository. We modified the provided code slightly to integrate it into our experiments.
  - `lsqrSOL`: We used the Matlab implementation of LSQR available at [https://web.stanford.edu/group/SOL/software/lsqr/](https://web.stanford.edu/group/SOL/software/lsqr/), and modified it to include the stopping criteria in equation (26) of the HIPNEX paper.
  - `minres_sol`: We implement our version of the MINRES algorithm, using the Matlab implementation available in the [Wikipedia page](https://en.wikipedia.org/wiki/Minimal_residual_method), as reference.