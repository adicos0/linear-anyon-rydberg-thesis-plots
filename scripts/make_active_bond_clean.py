"""
make_active_bond_clean.py
---------------------------------------------------------------------
Builds figures/rydberg_active_bond_clean.pdf from figures/rydberg_active_bond.pdf.

The source schematic is ~907 pt wide and carries a layer of 5 pt annotations
that are illegible at any scale the thesis page allows.  Per the supervisor's
comment ("eliminate/shorten the text that was too small to be seen; the caption
is already very detailed") this script strips them and keeps the drawing:

  * every text span below 7 pt in the diagram area,
  * the 14.3 pt composite title (the LaTeX caption repeats it),
  * the explanatory boxes at the foot of the source page (y > 468 pt),
  * the now-empty frame of the "Sign of phi" note box,

then crops to what is left.  Line art (levels, arrows, dotted guides) is kept
throughout, so the three panels are unchanged as drawings.
"""
import fitz

SRC = 'figures/rydberg_active_bond.pdf'
DST = 'figures/rydberg_active_bond_clean.pdf'
SIGN_BOX = fitz.Rect(812, 72, 902, 153)      # frame of the removed note box
CAPTION_Y = 468.0                            # foot-of-page explanatory boxes

d = fitz.open(SRC)
p = d[0]
W, H = p.rect.x1, p.rect.y1

# ---- stage 1: text only, line art preserved --------------------------------
# Filtering is done LINE by line, not span by span: the phase labels of panel (B)
# are a 7 pt base with a 5 pt superscript, so a span-level cut would strip the
# exponent and leave a meaningless "e".  A line goes if any part of it is below
# 6 pt.  The 6 pt subscripts of E_0 +- Delta_n therefore survive with their base.
n = 0
for b in p.get_text('dict')['blocks']:
    for l in b.get('lines', []):
        szs = [sp['size'] for sp in l['spans']]
        r = fitz.Rect(l['bbox'])
        if r.y1 > CAPTION_Y or min(szs) < 5.5 or max(szs) > 12:
            p.add_redact_annot(r); n += 1
p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                   graphics=fitz.PDF_REDACT_LINE_ART_NONE)
print('removed %d text lines' % n)

# ---- stage 2: the orphaned note-box frame, line art included ---------------
p.add_redact_annot(SIGN_BOX)
p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                   graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED)

# ---- crop to remaining content ---------------------------------------------
items = []
for dr in p.get_drawings():
    r = dr['rect']
    if r.y1 <= CAPTION_Y + 4 and r.width < W * 0.98 and r.height < H * 0.98:
        items.append(r)
for b in p.get_text('dict')['blocks']:
    for l in b.get('lines', []):
        for sp in l['spans']:
            r = fitz.Rect(sp['bbox'])
            if r.y1 <= CAPTION_Y + 4:
                items.append(r)
m = 8.0
box = fitz.Rect(max(0, min(r.x0 for r in items) - m), max(0, min(r.y0 for r in items) - m),
                min(W, max(r.x1 for r in items) + m), min(H, max(r.y1 for r in items) + m))
p.set_cropbox(box)
d.save(DST)

c = fitz.open(DST)[0]
sizes = sorted({round(sp['size'], 1) for b in c.get_text('dict')['blocks']
                for l in b.get('lines', []) for sp in l['spans']})
print('%s: %.1f x %.1f pt, surviving font sizes %s' % (DST, c.rect.width, c.rect.height, sizes))

# ---- also emit the three panels as separate crops -------------------------
# Boxes are in SOURCE page coordinates, measured from the redacted drawing.
# Panel C starts at 627 rather than at its content edge (623.6) so that the
# right border of panel B's header box (625.3) is excluded; the 3.4 pt of
# dotted guide lost that way is invisible.
PANELS = {'A': (41.8,  40.1, 362.4, 406.4),
          'B': (396.7, 40.1, 631.3, 450.6),
          'C': (627.0,  0.4, 900.6, 458.3)}
for name, (bx0, by0, bx1, by1) in PANELS.items():
    q = fitz.open(SRC)
    pg = q[0]
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            szs = [sp['size'] for sp in l['spans']]
            r = fitz.Rect(l['bbox'])
            if r.y1 > CAPTION_Y or min(szs) < 5.5 or max(szs) > 12:
                pg.add_redact_annot(r)
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                        graphics=fitz.PDF_REDACT_LINE_ART_NONE)
    pg.add_redact_annot(SIGN_BOX)
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                        graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED)
    pg.set_cropbox(fitz.Rect(bx0, by0, bx1, by1))
    out = 'figures/rydberg_active_bond_%s.pdf' % name
    q.save(out)
    print('%s: %.0f x %.0f pt' % (out, bx1 - bx0, by1 - by0))
