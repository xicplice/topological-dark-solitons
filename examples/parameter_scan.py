"""
Parameter scan over (lambda_eps, eps_inf) for the soliton model.

Usage (from repository root):

    python -m examples.parameter_scan --help
"""

from __future__ import annotations

import argparse

import numpy as np

from src.helmholtz_solver import SolitonParameters, defect_mass


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Scan soliton parameter space.")
    p.add_argument("--R", type=float, default=1.0, help="Macroscopic scale R.")
    p.add_argument("--n1", type=int, default=12, help="Topological integer n1.")
    p.add_argument("--n2", type=int, default=93, help="Topological integer n2.")

    p.add_argument(
        "--eps-inf-min",
        type=float,
        default=1e-7,
        help="Minimum epsilon_infty.",
    )
    p.add_argument(
        "--eps-inf-max",
        type=float,
        default=1e-5,
        help="Maximum epsilon_infty.",
    )
    p.add_argument(
        "--n-eps",
        type=int,
        default=5,
        help="Number of eps_inf samples.",
    )

    p.add_argument(
        "--lambda-factor-min",
        type=float,
        default=0.5,
        help="Minimum factor multiplying R / sqrt(n1 n2).",
    )
    p.add_argument(
        "--lambda-factor-max",
        type=float,
        default=2.0,
        help="Maximum factor multiplying R / sqrt(n1 n2).",
    )
    p.add_argument(
        "--n-lambda",
        type=int,
        default=5,
        help="Number of lambda_eps samples.",
    )

    p.add_argument(
        "--r-core",
        type=float,
        default=0.05,
        help="Core radius r_core.",
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()

    base_lambda = args.R / np.sqrt(args.n1 * args.n2)
    lambda_factors = np.linspace(
        args.lambda_factor_min, args.lambda_factor_max, args.n_lambda
    )
    eps_values = np.linspace(args.eps_inf_min, args.eps_inf_max, args.n_eps)

    print("# R = {:.3g}, n1 = {}, n2 = {}, base_lambda = {:.3e}".format(
        args.R, args.n1, args.n2, base_lambda
    ))
    print("# r_core = {:.3g}".format(args.r_core))
    print("# Columns: lambda_factor, lambda_eps, eps_inf, defect_mass")

    for lf in lambda_factors:
        lam = lf * base_lambda
        for eps_inf in eps_values:
            p = SolitonParameters(
                R=args.R,
                n1=args.n1,
                n2=args.n2,
                eps_inf=eps_inf,
                lambda_eps=lam,
                r_core=args.r_core,
                use_SI=False,
            )
            dm = defect_mass(p)
            print(
                "{:12.5g} {:12.5e} {:12.5e} {:12.5e}".format(
                    lf, lam, eps_inf, dm
                )
            )


if __name__ == "__main__":
    main()
