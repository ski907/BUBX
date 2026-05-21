import io
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Image as RLImage, HRFlowable,
)
from reportlab.lib.enums import TA_CENTER

from geometry.parse_results import parse_results
from physics import air
from conversions import convert


def _chart_orifice_flows(solved_geom, units='English'):
    results = parse_results(solved_geom)
    if units == 'SI':
        flows = results.orifice_flows_SCMM()
        offsets = results.orifice_positions()
        ylabel = 'Orifice Airflow (SCMM @ 1atm, 20°C)'
        xlabel = 'Orifice Offset (m)'
    else:
        flows = results.orifice_flows_SCFM()
        offsets = [convert.m_to_ft(o) for o in results.orifice_positions()]
        ylabel = 'Orifice Airflow (SCFM @ 1atm, 68°F)'
        xlabel = 'Orifice Offset (ft)'

    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.scatter(offsets, flows, color='steelblue', s=60, zorder=3)
    ax.axhline(np.mean(flows), color='gray', linestyle='--', linewidth=1,
               label=f'Mean: {np.mean(flows):.4f}')
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title('Orifice Flow Distribution', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return buf


def _chart_surface_velocities(solved_geom, units='English'):
    results = parse_results(solved_geom)
    if units == 'SI':
        vels = results.horizontal_surface_vel_haehnel2016_at_orifices()
        offsets = results.orifice_positions()
        ylabel = 'Surface Velocity (m/s)'
        xlabel = 'Orifice Offset (m)'
    else:
        vels = [convert.m_to_ft(v) for v in results.horizontal_surface_vel_haehnel2016_at_orifices()]
        offsets = [convert.m_to_ft(o) for o in results.orifice_positions()]
        ylabel = 'Surface Velocity (ft/s)'
        xlabel = 'Orifice Offset (ft)'

    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.scatter(offsets, vels, color='tomato', s=60, zorder=3)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title('Horizontal Surface Velocity Profile (Haehnel 2016)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return buf


def _base_table_style(header_bg=colors.HexColor('#1f4e79')):
    return TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), header_bg),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 10),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 9),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f0f5ff')]),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#c0c0c0')),
        ('TOPPADDING', (0, 0), (-1, -1), 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
        ('LEFTPADDING', (0, 0), (-1, -1), 8),
        ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
    ])


def _chart_diffuser_schematic(inputs: dict) -> io.BytesIO:
    """
    Draw a dimensioned side-view schematic of the diffuser system.
    Coordinate system: x = offset along diffuser (ft), y = elevation
    where y=0 is the diffuser pipe and y=water_depth_ft is the surface.
    """
    seg_len    = inputs.get('segment_length_ft', 100.0)
    n_or       = int(inputs.get('number_of_orifices', 10))
    or_dia     = inputs.get('orifice_diameter_in', 0.625)
    pipe_dia   = inputs.get('pipe_diameter_in', 3.0)
    depth      = inputs.get('water_depth_ft', 30.0)
    sup_len    = inputs.get('supply_pipe_length_ft', 100.0)
    sup_dia    = inputs.get('supply_pipe_diameter_in', 3.0)
    spacing    = seg_len / (n_or - 1) if n_or > 1 else seg_len

    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    ax.set_aspect('auto')
    ax.axis('off')

    # ── Water body ───────────────────────────────────────────────────────────
    x_left  = -seg_len * 0.06
    x_right =  seg_len * 1.06
    ax.fill_between([x_left, x_right], [0, 0], [depth, depth],
                    color='#cce5f6', alpha=0.45, zorder=0)
    ax.plot([x_left, x_right], [depth, depth],
            color='#2196F3', linewidth=1.8, zorder=3)
    ax.text(seg_len / 2, depth + depth * 0.03,
            'Water Surface', ha='center', va='bottom',
            fontsize=8, color='#1565C0', style='italic')

    # ── Seabed / bottom hatching ─────────────────────────────────────────────
    ax.fill_between([x_left, x_right], [-depth * 0.06, -depth * 0.06], [0, 0],
                    color='#b8a98a', alpha=0.5, zorder=0)
    ax.plot([x_left, x_right], [0, 0], color='#6d5a3e', linewidth=1, zorder=1)

    # ── Supply line (vertical, entering diffuser at left end) ─────────────────
    sup_show = min(depth * 0.35, sup_len * 0.12)   # portion visible in frame
    ax.plot([0, 0], [0, depth * 0.88],
            color='#444444', linewidth=4.5, solid_capstyle='round', zorder=4)
    ax.annotate('', xy=(0, depth * 0.88), xytext=(0, depth * 0.95),
                arrowprops=dict(arrowstyle='->', color='#444444', lw=1.5))
    ax.text(-seg_len * 0.025, depth * 0.55,
            f'Supply line\n{sup_dia:.2f}" ø\n{sup_len:.0f} ft',
            ha='right', va='center', fontsize=7.5, color='#333333',
            linespacing=1.4)

    # ── Diffuser pipe ────────────────────────────────────────────────────────
    ax.plot([0, seg_len], [0, 0],
            color='#222222', linewidth=5.5, solid_capstyle='butt', zorder=4)
    # Dead-end plug
    ax.plot([seg_len, seg_len], [-depth * 0.015, depth * 0.015],
            color='#222222', linewidth=5, solid_capstyle='butt', zorder=5)

    # ── Orifices ─────────────────────────────────────────────────────────────
    or_h   = depth * 0.10                  # nozzle stub height
    bub_h  = depth * 0.72                  # max bubble height
    x_ors  = [i * spacing for i in range(n_or)]

    # Decide how many orifice labels to show (avoid clutter)
    label_every = max(1, n_or // 6)

    for idx, xo in enumerate(x_ors):
        # Orifice marker directly on the pipe
        ax.plot(xo, 0, 'o',
                color='#e65100', markersize=5, zorder=5)

        # Bubble trail (3 dots rising from the pipe)
        # First orifice bubbles drift right to clear the supply line
        bubble_x_offset = spacing * 0.15 if idx == 0 else 0.0
        for frac in [0.25, 0.50, 0.78]:
            by = frac * bub_h
            bx = xo + bubble_x_offset + (frac - 0.5) * spacing * 0.08
            r  = max(1.5, 6 * frac)
            ax.plot(bx, by, 'o', color='#90caf9',
                    markersize=r, alpha=0.55, zorder=2)

    # ── Dimension: water depth (right side) ──────────────────────────────────
    dx = seg_len * 1.09
    ax.annotate('', xy=(dx, 0), xytext=(dx, depth),
                arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.1))
    ax.text(dx + seg_len * 0.012, depth / 2,
            f'{depth:.1f} ft\nwater depth',
            ha='left', va='center', fontsize=7.5, color='#222222',
            linespacing=1.4)

    # ── Dimension: total diffuser length (below pipe) ─────────────────────────
    dy1 = -depth * 0.10
    ax.plot([0, 0], [0, dy1], color='#666666', linewidth=0.7, linestyle='--')
    ax.plot([seg_len, seg_len], [0, dy1],
            color='#666666', linewidth=0.7, linestyle='--')
    ax.annotate('', xy=(0, dy1), xytext=(seg_len, dy1),
                arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.1))
    ax.text(seg_len / 2, dy1 - depth * 0.025,
            f'Diffuser length = {seg_len:.1f} ft   |   '
            f'Pipe {pipe_dia:.2f}" ø',
            ha='center', va='top', fontsize=7.5, color='#222222')

    # ── Dimension: orifice spacing (between first two) ───────────────────────
    if n_or > 1:
        dy2 = -depth * 0.22
        x0, x1 = x_ors[0], x_ors[1]
        ax.plot([x0, x0], [dy1, dy2], color='#888888', linewidth=0.7, linestyle='--')
        ax.plot([x1, x1], [dy1, dy2], color='#888888', linewidth=0.7, linestyle='--')
        ax.annotate('', xy=(x0, dy2), xytext=(x1, dy2),
                    arrowprops=dict(arrowstyle='<->', color='#555555', lw=1.0))
        ax.text((x0 + x1) / 2, dy2 - depth * 0.025,
                f'Spacing = {spacing:.1f} ft   |   '
                f'Orifice {or_dia:.4f}" ø  x{n_or}',
                ha='center', va='top', fontsize=7.5, color='#444444')

    # ── Axes limits and title ─────────────────────────────────────────────────
    y_bot = -depth * 0.35
    y_top =  depth * 1.12
    ax.set_xlim(x_left - seg_len * 0.12, x_right + seg_len * 0.20)
    ax.set_ylim(y_bot, y_top)
    ax.set_title('Diffuser System Schematic  (not to scale)',
                 fontsize=10, fontweight='bold', pad=4)

    plt.tight_layout(pad=0.4)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=160, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return buf


def generate_pdf_report(solved_geom, inputs: dict, air_pressure: float,
                        water_pressure: float, air_temp: float) -> bytes:
    """
    Generate a summary PDF. Returns raw bytes suitable for st.download_button.

    inputs keys:
        pipe_type, pipe_diameter_in, segment_length_ft, number_of_orifices,
        orifice_diameter_in, supply_pipe_type, supply_pipe_diameter_in,
        supply_pipe_length_ft, bc_method, water_depth_ft
    """
    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf, pagesize=letter,
        leftMargin=0.75 * inch, rightMargin=0.75 * inch,
        topMargin=0.75 * inch, bottomMargin=0.75 * inch,
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'ReportTitle', parent=styles['Title'],
        fontSize=20, textColor=colors.HexColor('#1f4e79'), spaceAfter=2,
    )
    meta_style = ParagraphStyle(
        'Meta', parent=styles['Normal'],
        fontSize=9, textColor=colors.HexColor('#666666'), spaceAfter=10,
    )
    section_style = ParagraphStyle(
        'Section', parent=styles['Heading2'],
        fontSize=12, textColor=colors.HexColor('#1f4e79'),
        spaceBefore=14, spaceAfter=6,
    )
    footer_style = ParagraphStyle(
        'Footer', parent=styles['Normal'],
        fontSize=8, textColor=colors.HexColor('#999999'), alignment=TA_CENTER,
    )

    story = []

    # Header
    run_name = inputs.get('run_name', '').strip()
    story.append(Paragraph('BUBX Air Demand Calculator', title_style))
    if run_name:
        run_name_style = ParagraphStyle(
            'RunName', parent=styles['Normal'],
            fontSize=13, textColor=colors.HexColor('#1f4e79'),
            spaceAfter=3, fontName='Helvetica-Bold',
        )
        story.append(Paragraph(run_name, run_name_style))
    story.append(Paragraph(
        f'Summary Report  |  Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}  |  v0.2',
        meta_style,
    ))
    story.append(HRFlowable(width='100%', thickness=2, color=colors.HexColor('#1f4e79')))
    story.append(Spacer(1, 10))

    # Input Parameters
    story.append(Paragraph('Input Parameters', section_style))

    gauge_psi = convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure))
    gauge_mpa = convert.pressure_to_gauge(air_pressure) / 1e6

    inputs_data = [
        ['Parameter', 'Value'],
        ['Diffuser Pipe Material', inputs.get('pipe_type', '-')],
        ['Diffuser Pipe Diameter', f"{inputs.get('pipe_diameter_in', 0):.2f} in"],
        ['Diffuser Segment Length', f"{inputs.get('segment_length_ft', 0):.1f} ft"],
        ['Number of Orifices', str(inputs.get('number_of_orifices', '-'))],
        ['Orifice Spacing', (
            f"{inputs['segment_length_ft'] / (inputs['number_of_orifices'] - 1):.2f} ft"
            if inputs.get('number_of_orifices', 1) > 1 else '-'
        )],
        ['Orifice Diameter', f"{inputs.get('orifice_diameter_in', 0):.4f} in"],
        ['Supply Pipe Material', inputs.get('supply_pipe_type', '-')],
        ['Supply Pipe Diameter', f"{inputs.get('supply_pipe_diameter_in', 0):.2f} in"],
        ['Supply Pipe Length', f"{inputs.get('supply_pipe_length_ft', 0):.1f} ft"],
        ['Boundary Condition', inputs.get('bc_method', '-')],
        ['Water Depth', f"{inputs.get('water_depth_ft', 0):.1f} ft"],
        ['Air Pressure (gauge)', f"{gauge_psi:.2f} psi  /  {gauge_mpa:.3f} MPa"],
        ['Air Temperature', f"{air_temp:.0f} °C"],
    ]

    t_inputs = Table(inputs_data, colWidths=[3.0 * inch, 3.5 * inch])
    t_inputs.setStyle(_base_table_style())
    story.append(t_inputs)
    story.append(Spacer(1, 10))

    # Diffuser schematic
    story.append(RLImage(_chart_diffuser_schematic(inputs),
                         width=6.5 * inch, height=2.75 * inch))
    story.append(Spacer(1, 12))

    # Flow Results
    story.append(Paragraph('Flow Rates & System Parameters', section_style))

    results = parse_results(solved_geom)
    mdot = results.total_mdot()
    atm_p = convert.pressure_to_absolute(0)
    rho_68f = air.rho_air(atm_p, T=convert.F_to_C(68))

    flow_si = convert.CMS_to_CMM(air.Q(mdot, rho_68f))
    flow_imp = convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot, rho_68f)))
    q_per_len = results.airflow_per_unit_length()
    vel_ms = results.horizontal_surface_vel_haehnel2016()
    vel_fts = convert.m_to_ft(vel_ms)

    flow_data = [
        ['Parameter', 'SI', 'Imperial'],
        ['Total Flow Rate', f'{flow_si:.2f} SCMM', f'{flow_imp:.2f} SCFM'],
        ['Air Pressure (gauge)', f'{gauge_mpa:.3f} MPa', f'{gauge_psi:.2f} psi'],
        ['Airflow Per Unit Length', f'{q_per_len:.3f} SCMM/m', f'{convert.CMM_to_CFM(q_per_len) * convert.ft_to_m(1):.3f} SCFM/ft'],
        ['Surface Velocity (Haehnel 2016)', f'{vel_ms:.3f} m/s', f'{vel_fts:.3f} ft/s'],
    ]

    t_flow = Table(flow_data, colWidths=[2.75 * inch, 1.875 * inch, 1.875 * inch])
    t_flow.setStyle(_base_table_style())
    story.append(t_flow)
    story.append(Spacer(1, 12))

    # Performance Metrics
    story.append(Paragraph('Performance Metrics', section_style))

    cu = results.coefficient_of_uniformity()
    odr = results.orifice_to_diffuser_area_ratio()
    cu_pass = cu >= 0.9
    odr_pass = odr <= 0.25

    perf_style = _base_table_style()
    for row_idx, passing in enumerate([cu_pass, odr_pass], start=1):
        bg = colors.HexColor('#d4edda') if passing else colors.HexColor('#fff3cd')
        tc = colors.HexColor('#155724') if passing else colors.HexColor('#856404')
        perf_style.add('BACKGROUND', (3, row_idx), (3, row_idx), bg)
        perf_style.add('TEXTCOLOR', (3, row_idx), (3, row_idx), tc)
        perf_style.add('FONTNAME', (3, row_idx), (3, row_idx), 'Helvetica-Bold')

    perf_data = [
        ['Metric', 'Value', 'Target', 'Status'],
        ['Coefficient of Uniformity', f'{cu:.3f}', '>= 0.900', 'PASS' if cu_pass else 'WARNING'],
        ['Orifice / Diffuser Area Ratio', f'{odr:.3f}', '<= 0.250', 'PASS' if odr_pass else 'WARNING'],
    ]

    t_perf = Table(perf_data, colWidths=[2.5 * inch, 1.25 * inch, 1.25 * inch, 1.5 * inch])
    t_perf.setStyle(perf_style)
    story.append(t_perf)
    story.append(Spacer(1, 16))

    # Charts
    story.append(HRFlowable(width='100%', thickness=1, color=colors.HexColor('#dddddd')))
    story.append(Paragraph('Orifice Flow Distribution', section_style))
    story.append(RLImage(_chart_orifice_flows(solved_geom), width=6.5 * inch, height=3.25 * inch))
    story.append(Spacer(1, 12))

    story.append(Paragraph('Horizontal Surface Velocity Profile', section_style))
    story.append(RLImage(_chart_surface_velocities(solved_geom), width=6.5 * inch, height=3.25 * inch))

    # Footer
    story.append(Spacer(1, 14))
    story.append(HRFlowable(width='100%', thickness=0.5, color=colors.HexColor('#dddddd')))
    story.append(Paragraph(
        'BUBX v0.2  |  Air Demand Calculator for Submerged Manifolds  |  '
        'Surface velocity per Haehnel (2016)',
        footer_style,
    ))

    doc.build(story)
    buf.seek(0)
    return buf.getvalue()
