"""
compare_models.py
-----------------
Compares the BUBX compressible-isentropic orifice model against the
BUB300 Bernoulli / Torricelli model across:

  1. Single-orifice mass-flow sweep over a range of pressure ratios
     (isolates pure orifice physics; no pipe solver involved)

  2. Full diffuser system for four representative operating cases

  3. Per-orifice flow distribution for the first test case

Run with:  python compare_models.py
Outputs:   console table  +  model_comparison.png
"""

import math
import copy
import io
import contextlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from geometry.orifice import Orifice
from geometry.segment import Segment
from geometry.segments import generate_equal_spaced_uniform, add_supply_line
from solver import downstream_solve
from geometry.parse_results import parse_results
from physics import air
from conversions import convert


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

PATM_BUB300    = 98100    # Pa  – BUB300 regional atmosphere (Hanover, NH)
DENAIR0_BUB300 = 1.294    # kg/m³ – BUB300 reference density @ 0 °C
GAMMA          = 1.4
RC             = (2 / (GAMMA + 1)) ** (GAMMA / (GAMMA - 1))  # ≈ 0.528 choke criterion
ROUGHNESS_SS   = 2.5e-5   # m – stainless steel absolute roughness


# ─────────────────────────────────────────────────────────────────────────────
# BUB300-style orifice
# ─────────────────────────────────────────────────────────────────────────────

class BUB300Orifice(Orifice):
    """
    Drop-in replacement for BUBX Orifice using the BUB300 Bernoulli model.

    BUB300.f lines 496-503:
        P_diff  = P_air_gauge - P_water_gauge          (differential)
        Vout    = sqrt(2 * P_diff / Deno)              (Torricelli)
        mdot    = Cdis * Aorf * Deno * Vout

    Exit density (Deno) options controlled by exp_type:
        0 = incompressible  – density = local pipe density
        1 = adiabatic       – Deno = Den * (P_water/P_air)^(1/γ)
        2 = isothermal      – Deno = Denair0 * P_water_abs / Patm  ← BUB300 default

    Inherits from Orifice so isinstance() checks in the solver still pass.
    """

    def __init__(self, diameter: float, Cdis: float = 1.0, exp_type: int = 2,
                 Patm: float = PATM_BUB300, Denair0: float = DENAIR0_BUB300):
        super().__init__(diameter)
        self.Cdis     = Cdis
        self.exp_type = exp_type
        self._Patm    = Patm
        self._Denair0 = Denair0

    @property
    def mdot(self) -> float:
        P_air   = self.upstream_pressure
        P_water = self.downstream_pressure
        P_diff  = P_air - P_water
        if P_diff <= 0.0:
            return 0.0

        # Local pipe density (isothermal ideal gas scaled from reference)
        Den_pipe = self._Denair0 * P_air / self._Patm

        if self.exp_type == 0:
            Deno = Den_pipe
        elif self.exp_type == 1:
            Deno = Den_pipe * (P_water / P_air) ** (1.0 / GAMMA)
        else:                                          # isothermal (exp_type=2)
            Deno = self._Denair0 * P_water / self._Patm

        Vout = math.sqrt(2.0 * P_diff / Deno)
        return self.Cdis * self.area * Deno * Vout


# ─────────────────────────────────────────────────────────────────────────────
# Model registry
# ─────────────────────────────────────────────────────────────────────────────

MODELS = [
    {
        'name':   'BUBX – isentropic, Cc=0.90',
        'cls':    Orifice,
        'kwargs': {'contraction_coeff': 0.90},
        'color':  '#2166AC',
        'ls':     '-',
        'marker': 'o',
    },
    {
        'name':   'BUB300 – isothermal, Cdis=1.00 (as-coded)',
        'cls':    BUB300Orifice,
        'kwargs': {'Cdis': 1.00, 'exp_type': 2},
        'color':  '#D6604D',
        'ls':     '--',
        'marker': 's',
    },
    {
        'name':   'BUB300 – isothermal, Cdis=0.63 (per comment)',
        'cls':    BUB300Orifice,
        'kwargs': {'Cdis': 0.63, 'exp_type': 2},
        'color':  '#F4A582',
        'ls':     '-.',
        'marker': '^',
    },
    {
        'name':   'BUB300 – isothermal, Cdis=0.90 (matched Cc)',
        'cls':    BUB300Orifice,
        'kwargs': {'Cdis': 0.90, 'exp_type': 2},
        'color':  '#4DAC26',
        'ls':     ':',
        'marker': 'D',
    },
]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def swap_orifices(geom: Segment, OrificeClass, **kwargs) -> Segment:
    """
    Return a new Segment with every Orifice instance replaced by
    OrificeClass(diameter, **kwargs), preserving elevation attributes.
    """
    new_features = []
    for feature in geom.features:
        if isinstance(feature, Orifice):
            new_o = OrificeClass(feature.diameter, **kwargs)
            new_o.elevation_above_datum = feature.elevation_above_datum
            new_o.datum = feature.datum
            new_features.append(new_o)
        else:
            new_features.append(copy.deepcopy(feature))
    return Segment(new_features)


def build_geometry(pipe_dia_in: float, seg_len_ft: float, n_orifices: int,
                   orifice_dia_in: float, supply_len_ft: float) -> Segment:
    return add_supply_line(
        generate_equal_spaced_uniform(
            convert.ft_to_m(seg_len_ft),
            ROUGHNESS_SS,
            convert.in_to_m(pipe_dia_in),
            n_orifices,
            convert.in_to_m(orifice_dia_in),
        ),
        convert.ft_to_m(supply_len_ft),
        ROUGHNESS_SS,
        convert.in_to_m(pipe_dia_in),
    )


@contextlib.contextmanager
def _quiet():
    """Suppress the print(m) bisection trace inside downstream_solve."""
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def run_solver(geom: Segment, P_air: float, P_water: float, air_temp: float = 0):
    with _quiet():
        solved = downstream_solve(geom, P_air, P_water, 0.001, air_temp)
    return parse_results(solved)


def _to_scfm(mdot_kgs: float, rho_std: float) -> float:
    return convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot_kgs, rho_std)))


# ─────────────────────────────────────────────────────────────────────────────
# 1. Single-orifice sweep
# ─────────────────────────────────────────────────────────────────────────────

def single_orifice_sweep(orifice_dia_in: float = 0.625,
                         water_depth_ft: float = 30.0,
                         n_points: int = 100,
                         max_ratio: float = 3.0):
    """
    Sweep P_air/P_water from 1.02 to max_ratio without the pipe solver.
    Returns (ratios array, {model_name: mdot_g_per_s list}, choke_ratio).
    """
    d_m     = convert.in_to_m(orifice_dia_in)
    P_water = convert.pressure_to_absolute(
                  convert.H_m_to_Pa(convert.ft_to_m(water_depth_ft)))
    ratios  = np.linspace(1.02, max_ratio, n_points)
    results = {m['name']: [] for m in MODELS}

    for ratio in ratios:
        P_air = P_water * ratio
        for m in MODELS:
            o = m['cls'](d_m, **m['kwargs'])
            o.upstream_pressure   = P_air
            o.downstream_pressure = P_water
            o.temp = 0
            results[m['name']].append(o.mdot * 1000)   # g/s

    choke_ratio = 1.0 / RC   # P_air/P_water at onset of choked flow
    return ratios, results, choke_ratio


# ─────────────────────────────────────────────────────────────────────────────
# 2. Full-system comparison
# ─────────────────────────────────────────────────────────────────────────────

# (label, pipe_dia_in, seg_len_ft, n_orifices, orifice_dia_in, supply_len_ft, depth_ft, press_psi)
TEST_CASES = [
    ('30 ft | 90 psi',   3.0, 100, 10, 0.625, 100, 30,  90),
    ('45 ft | 90 psi',   3.0, 100, 10, 0.625, 100, 45,  90),
    ('200 ft diffuser',  3.0, 200, 20, 0.625, 100, 30,  90),
    ('30 ft | 120 psi',  3.0, 100, 10, 0.625, 100, 30, 120),
]


def run_full_system():
    """
    Run MODELS × TEST_CASES.
    Returns (summary DataFrame, orifice-flow distribution dict for case 0).
    """
    atm_p   = convert.pressure_to_absolute(0)
    rho_std = air.rho_air(atm_p, T=convert.F_to_C(68))

    summary_rows = []
    distribution = {}   # {model_name: [SCFM per orifice]} for TEST_CASES[0]

    for case_idx, (label, pipe_dia_in, seg_len_ft, n_orifices, orifice_dia_in,
                   supply_len_ft, depth_ft, press_psi) in enumerate(TEST_CASES):

        base_geom = build_geometry(pipe_dia_in, seg_len_ft, n_orifices,
                                   orifice_dia_in, supply_len_ft)
        P_air   = convert.pressure_to_absolute(convert.psi_to_Pa(press_psi))
        P_water = convert.pressure_to_absolute(
                      convert.H_m_to_Pa(convert.ft_to_m(depth_ft)))

        print(f'  Case: {label}')

        for m in MODELS:
            geom = swap_orifices(base_geom, m['cls'], **m['kwargs'])
            try:
                r          = run_solver(geom, P_air, P_water)
                total_scfm = _to_scfm(r.total_mdot(), rho_std)
                cu         = r.coefficient_of_uniformity()
                odr        = r.orifice_to_diffuser_area_ratio()

                summary_rows.append({
                    'Case':             label,
                    'Model':            m['name'],
                    'Total (SCFM)':     round(total_scfm, 1),
                    'CU':               round(cu, 3),
                    'Aorf/Apipe':       round(odr, 3),
                })

                if case_idx == 0:
                    distribution[m['name']] = r.orifice_flows_SCFM()

            except Exception as exc:
                summary_rows.append({
                    'Case':         label,
                    'Model':        m['name'],
                    'Total (SCFM)': f'ERR: {exc}',
                    'CU':           '-',
                    'Aorf/Apipe':   '-',
                })

    return pd.DataFrame(summary_rows), distribution


# ─────────────────────────────────────────────────────────────────────────────
# 3. Plots
# ─────────────────────────────────────────────────────────────────────────────

def plot_all(sweep_ratios, sweep_results, choke_ratio,
             summary_df, distribution, orifice_dia_in=0.625):

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    fig.subplots_adjust(wspace=0.36, left=0.06, right=0.97, top=0.88, bottom=0.14)

    # ── Panel A: single-orifice sweep ────────────────────────────────────────
    ax = axes[0]
    for m in MODELS:
        ax.plot(sweep_ratios, sweep_results[m['name']],
                color=m['color'], linestyle=m['ls'], linewidth=1.8,
                label=m['name'])
    ax.axvline(choke_ratio, color='gray', linestyle=':', linewidth=1.2,
               label=f'BUBX choke onset  (r = {choke_ratio:.2f})')
    ax.set_xlabel('P$_{air}$ / P$_{water}$  (absolute pressures)', fontsize=10)
    ax.set_ylabel('Orifice mass flow  (g/s)', fontsize=10)
    ax.set_title(f'A.  Single-orifice sweep\n'
                 f'{orifice_dia_in}" ø orifice, 30 ft water depth', fontsize=10)
    ax.legend(fontsize=7.5, loc='upper left')
    ax.grid(True, alpha=0.3)

    # ── Panel B: total system flow – grouped bar chart ────────────────────────
    ax = axes[1]
    n_cases  = len(TEST_CASES)
    n_models = len(MODELS)
    bar_w    = 0.18
    x_pos    = np.arange(n_cases)

    for j, m in enumerate(MODELS):
        flows = []
        for label, *_ in TEST_CASES:
            mask = (summary_df['Case'] == label) & (summary_df['Model'] == m['name'])
            row  = summary_df[mask]
            if not row.empty:
                try:
                    flows.append(float(row.iloc[0]['Total (SCFM)']))
                except (ValueError, TypeError):
                    flows.append(0.0)
            else:
                flows.append(0.0)
        offset = (j - (n_models - 1) / 2.0) * bar_w
        ax.bar(x_pos + offset, flows, bar_w,
               color=m['color'], alpha=0.85, label=m['name'])

    ax.set_xticks(x_pos)
    ax.set_xticklabels([tc[0] for tc in TEST_CASES], fontsize=8.5)
    ax.set_ylabel('Total system flow  (SCFM @ 1 atm, 68 °F)', fontsize=10)
    ax.set_title('B.  Full-system total flow', fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3, axis='y')

    # ── Panel C: per-orifice distribution for TEST_CASES[0] ─────────────────
    ax = axes[2]
    if distribution:
        n_or    = len(next(iter(distribution.values())))
        x_or    = np.arange(1, n_or + 1)
        for m in MODELS:
            if m['name'] in distribution:
                ax.plot(x_or, distribution[m['name']],
                        color=m['color'], linestyle=m['ls'],
                        marker=m['marker'], markersize=5,
                        linewidth=1.6, label=m['name'])
    ax.set_xlabel('Orifice number  (upstream → downstream)', fontsize=10)
    ax.set_ylabel('Orifice flow  (SCFM)', fontsize=10)
    ax.set_title(f'C.  Per-orifice distribution\n'
                 f'({TEST_CASES[0][0]})', fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)

    fig.suptitle('BUBX vs BUB300 Orifice Model Comparison', fontsize=13)
    plt.savefig('model_comparison.png', dpi=150, bbox_inches='tight')
    print('\nPlot saved -> model_comparison.png')
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':

    SEP = '=' * 68

    print(SEP)
    print('  BUBX vs BUB300 Physics Comparison')
    print(SEP)

    # ── 1. Single-orifice ────────────────────────────────────────────────────
    print('\n[1/2]  Single-orifice sweep  (5/8" ø, 30 ft depth) ...')
    ratios, sweep_results, choke_ratio = single_orifice_sweep()
    print(f'       BUBX choked-flow onset at P_air/P_water = {choke_ratio:.3f}')

    # Spot-check table at a handful of pressure ratios
    spot = [1.1, 1.5, 2.0, round(choke_ratio, 2), 3.0]
    idx  = [int(np.argmin(np.abs(ratios - r))) for r in spot]
    spot_rows = []
    for i in idx:
        spot_rows.append(
            {'P_air/P_water': f'{ratios[i]:.2f}'}
            | {m['name']: f'{sweep_results[m["name"]][i]:.3f}' for m in MODELS}
        )
    print('\n  Mass flow per orifice (g/s) at selected pressure ratios:')
    print(pd.DataFrame(spot_rows).to_string(index=False))

    # ── 2. Full system ────────────────────────────────────────────────────────
    print(f'\n[2/2]  Full-system runs ({len(TEST_CASES)} cases × {len(MODELS)} models) ...')
    summary_df, distribution = run_full_system()

    print('\n  Summary:')
    pd.set_option('display.max_colwidth', 55)
    pd.set_option('display.width', 220)
    print(summary_df.to_string(index=False))

    # ── 3. Plots ──────────────────────────────────────────────────────────────
    print('\n  Generating plots ...')
    plot_all(ratios, sweep_results, choke_ratio, summary_df, distribution)

    print('\nDone.')
