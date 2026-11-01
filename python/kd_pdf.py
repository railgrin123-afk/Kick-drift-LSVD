"""
kd_pdf -- build the short PDF report from kd_report's tables and figures.

    python3 kd_pdf.py [--report reports]

Layout: summary, tables (checks first; the Distribution of S section opens
with the 4C table), figures, then method and caveats.
"""
import argparse, json, os
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                PageBreak, Image, KeepTogether)

FONT_DIR = '/usr/share/fonts/truetype/dejavu/'
pdfmetrics.registerFont(TTFont('DV', FONT_DIR + 'DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DVB', FONT_DIR + 'DejaVuSans-Bold.ttf'))
from reportlab.lib.fonts import addMapping
addMapping('DV', 0, 0, 'DV'); addMapping('DV', 1, 0, 'DVB'); addMapping('DV', 0, 1, 'DV'); addMapping('DV', 1, 1, 'DVB')

ss = getSampleStyleSheet()
H1 = ParagraphStyle('H1', parent=ss['Heading1'], fontName='DVB', fontSize=14, spaceAfter=6, spaceBefore=4)
H2 = ParagraphStyle('H2', parent=ss['Heading2'], fontName='DVB', fontSize=11, spaceAfter=4, spaceBefore=8)
B = ParagraphStyle('B', parent=ss['BodyText'], fontName='DV', fontSize=8.6, leading=11.4, spaceAfter=4)
SM = ParagraphStyle('SM', parent=B, fontSize=7.4, leading=9.4, textColor=colors.HexColor('#52514e'))
CELL = ParagraphStyle('CELL', parent=B, fontSize=6.6, leading=8, spaceAfter=0)
CELLH = ParagraphStyle('CELLH', parent=CELL, fontName='DVB')
TITLE = ParagraphStyle('T', parent=H1, fontSize=18, spaceAfter=2)


def table(hdr, rows, widths=None, pass_col=None):
    data = [[Paragraph(str(h), CELLH) for h in hdr]] + [[Paragraph(str(c), CELL) for c in r] for r in rows]
    t = Table(data, colWidths=widths, repeatRows=1)
    st = [('LINEBELOW', (0, 0), (-1, 0), 0.6, colors.HexColor('#52514e')),
          ('VALIGN', (0, 0), (-1, -1), 'TOP'),
          ('TOPPADDING', (0, 0), (-1, -1), 1.5), ('BOTTOMPADDING', (0, 0), (-1, -1), 1.5),
          ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f4f3f0')])]
    t.setStyle(TableStyle(st))
    return t


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument('--report', default=os.path.join(os.path.dirname(here), 'reports'))
    a = ap.parse_args()
    J = json.load(open(os.path.join(a.report, 'kd_tables.json')))
    F = json.load(open(os.path.join(a.report, 'kd_findings.json')))
    T, checks, figs, nseed, C = J['T'], J['checks'], J['figs'], J['nseed'], J['cycles']
    W = landscape(A4)[0] - 2.4 * cm
    doc = SimpleDocTemplate(os.path.join(a.report, 'kd_report.pdf'), pagesize=landscape(A4),
                            leftMargin=1.2 * cm, rightMargin=1.2 * cm, topMargin=1.0 * cm, bottomMargin=1.0 * cm,
                            title='Kick-drift LSVD: manifolds, crossing ages and the distribution of S',
                            author='Kick-drift LSVD analysis')
    S = []
    S.append(Paragraph('Kick-drift LSVD — manifolds, crossing ages and the distribution of S', TITLE))
    S.append(Paragraph('%s; %d cycles; Python engine (python/kd_engine.py), multi-seed, mean ± sd over seeds. '
                       'Seed 0 reproduces the MATLAB rng(0) run exactly.' %
                       (', '.join('%s %d seeds' % (m, n) for m, n in nseed.items()), C), SM))
    S.append(Paragraph('Summary', H2))
    for f in F['summary']:
        S.append(Paragraph('• ' + f, B))

    # ---------------- tables
    S.append(Paragraph('Checks', H2))
    rows = [['<font color="%s"><b>%s</b></font>' % ('#008300' if ok else '#e34948', 'PASS' if ok else 'FAIL'), n, d]
            for n, ok, d in checks]
    S.append(table(['', 'check', 'detail'], rows, [1.3 * cm, 15 * cm, W - 16.3 * cm]))

    def tsec(key, title, note=None):
        hdr, rows = T[key]
        blk = [Paragraph(title, H2)]
        if note:
            blk.append(Paragraph(note, SM))
        blk.append(table(hdr, rows))
        S.append(KeepTogether(blk) if len(rows) < 25 else blk[0])
        if len(rows) >= 25:
            S.extend(blk[1:])

    S.append(PageBreak())
    S.append(Paragraph('Part 1 — launched lines of rays', H1))
    tsec('caustics', '1A. Caustic growth vs Kaplan Eq. (6) prediction α + β/4', F['notes'].get('caustics'))
    tsec('corr', '1A(iv). Correlation across x of the ray density ρ', F['notes'].get('corr'))
    tsec('1C', '1C. Do low-stretch rays make the branches? (launched lines, mean of 3 rows)', F['notes'].get('1C'))
    tsec('1C_lat', '1C control: the tracer lattice (area-filling, so its projection has no branches)')
    S.append(PageBreak())
    S.append(Paragraph('Part 2 — crossing ages', H1))
    tsec('pooled', '2A. Pooled age-bin hazards', F['notes'].get('pooled'))
    tsec('KM', '2B. Kaplan–Meier F(k) and the naive fraction', F['notes'].get('KM'))
    tsec('depth', '2C. Does "returned before" raise the re-entry hazard at fixed depth?', F['notes'].get('depth'))
    S.append(PageBreak())
    S.append(Paragraph('Distribution of S', H1))
    tsec('4C', '4C. Kaplan (PRL 89, 184103) log-normal tail test, fit N = 20..%d' % C, F['notes'].get('4C'))
    tsec('4C_seed0', '4C reproduction, seed 0')
    tsec('Ple', '4A. P(S ≤ s) (mean over seeds)')
    tsec('mom', '4A. Moments and quantiles of ln S')
    tsec('iid', '4D. i.i.d. random-matrix model vs the 2-D annealed run', F['notes'].get('iid'))
    tsec('iid_P', '4D. P(S ≤ s, N): i.i.d. vs 2-D annealed (z = difference in standard errors of the seed mean)')
    tsec('Lq', '4D. Generalised Lyapunov exponents L(q)', F['notes'].get('Lq'))
    tsec('4E', '4E. Robustness: threshold and metric', F['notes'].get('4E'))
    tsec('4F', '4F. Rarefaction r vs ln S on the lines', F['notes'].get('4F'))

    # ---------------- figures
    S.append(PageBreak())
    S.append(Paragraph('Figures', H1))
    from reportlab.lib.utils import ImageReader
    maxw, maxh = W, landscape(A4)[1] - 4.2 * cm
    for p, cap in figs:
        iw, ih = ImageReader(p).getSize()
        sc = min(maxw / iw, maxh / ih)
        S.append(KeepTogether([Image(p, iw * sc, ih * sc), Paragraph(cap, SM), Spacer(1, 6)]))

    # ---------------- method and caveats
    S.append(PageBreak())
    S.append(Paragraph('Method', H1))
    for f in F['method']:
        S.append(Paragraph(f, B))
    S.append(Paragraph('Caveats', H1))
    for f in F['caveats']:
        S.append(Paragraph('• ' + f, B))
    doc.build(S)
    print('written', os.path.join(a.report, 'kd_report.pdf'))


if __name__ == '__main__':
    main()
