#!/usr/bin/env python
"""Brand assets from the team logo (this is what `make brand` runs).

Input   branding/ratpac-logo-1600x1800.png   1600 x 1800 RGBA, opaque cream (#f4efe6) background, the round navy
                                             badge in the upper part and the wordmark "The Rat PAC" below
        site/data/headline.json["question"]  the question printed on the social-preview card
Output  site/assets/brand/badge-{16,32,64,192,512}.png   the badge disc, transparent outside the circle
        site/assets/brand/badge-180.png                  the disc on an opaque cream square (apple-touch-icon)
        site/favicon.ico                                 16 / 32 / 48
        site/assets/brand/ratpac-logo.png                the full logo at 800 x 900
        site/assets/brand/badge-og.png                   1200 x 630 social preview: badge, question, team, URL

The badge is detected, not assumed: within the top 1,250 rows, pixels whose RGB distance from the cream exceeds 40
are labelled (scipy.ndimage.label); the largest connected component is the badge ring, and its bounding box gives
the centre and radius, which must lie within 30 px of the expected centre (800, 690) and radius 490.

Usage: python tools/make_logo_assets.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "branding" / "ratpac-logo-1600x1800.png"
HEADLINE = ROOT / "site" / "data" / "headline.json"
BRAND = ROOT / "site" / "assets" / "brand"
FAVICON = ROOT / "site" / "favicon.ico"

CREAM = (0xF4, 0xEF, 0xE6)
NAVY = (0x1F, 0x2A, 0x3C)
ORANGE = (0xE8, 0x5A, 0x2B)
INK2 = (0x52, 0x51, 0x4E)
EXPECTED = {"cx": 800.0, "cy": 690.0, "r": 490.0}
TOLERANCE = 30          # px, on each of cx, cy and r
TOP_ROWS = 1250         # the badge lives above the wordmark
CREAM_DISTANCE = 40     # RGB distance beyond which a pixel is "not background"
MARGIN = 0.035          # crop margin around the disc, as a fraction of the radius (3-4 %)
MASK_RADIUS = 0.985     # circular alpha-mask radius, as a fraction of the crop half-size
FEATHER = 1.5           # px, linear alpha ramp at the mask edge
SIZE_WARN = 400_000     # bytes
FONT_DIR = Path("/usr/share/fonts/truetype/dejavu")
LANCZOS = Image.Resampling.LANCZOS
SITE_URL = "corygardner.github.io/motrpac_challenge"
TEAM_LINE = "The Rat PAC · MoTrPAC Hackathon 2026"


def detect_badge(im: Image.Image) -> dict:
    """Centre, radius and bounding box of the badge ring: the largest non-cream connected component in the top rows."""
    rgb = np.asarray(im.convert("RGB"), dtype=np.int32)[:TOP_ROWS]
    dist = np.sqrt(((rgb - np.array(CREAM, dtype=np.int32)) ** 2).sum(axis=2))
    mask = dist > CREAM_DISTANCE
    labels, n = ndimage.label(mask)
    if n == 0:
        raise SystemExit(f"no pixel in the top {TOP_ROWS} rows is farther than {CREAM_DISTANCE} from the cream")
    sizes = ndimage.sum(mask, labels, index=np.arange(1, n + 1))
    k = int(np.argmax(sizes)) + 1
    ys, xs = np.nonzero(labels == k)
    x0, x1, y0, y1 = int(xs.min()), int(xs.max()), int(ys.min()), int(ys.max())
    w, h = x1 - x0 + 1, y1 - y0 + 1
    return {"cx": (x0 + x1 + 1) / 2, "cy": (y0 + y1 + 1) / 2, "r": (w + h) / 4, "bbox": (x0, y0, x1, y1),
            "w": w, "h": h, "components": int(n), "pixels": int(sizes[k - 1])}


def crop_disc(im: Image.Image, g: dict) -> tuple[Image.Image, tuple[int, int, int, int]]:
    """Square crop around the disc with a small margin, alpha-masked to a feathered circle."""
    half = g["r"] * (1 + MARGIN)
    side = int(round(2 * half))
    x0, y0 = int(round(g["cx"] - half)), int(round(g["cy"] - half))
    box = (x0, y0, x0 + side, y0 + side)
    if x0 < 0 or y0 < 0 or box[2] > im.size[0] or box[3] > im.size[1]:
        raise SystemExit(f"crop box {box} leaves the {im.size} image")
    sq = np.array(im.crop(box))  # RGBA copy
    yy, xx = np.mgrid[:side, :side]
    d = np.hypot(xx + 0.5 - side / 2, yy + 0.5 - side / 2)
    alpha = np.clip((MASK_RADIUS * side / 2 - d) / FEATHER + 0.5, 0.0, 1.0)
    sq[..., 3] = np.round(alpha * 255).astype(np.uint8)
    return Image.fromarray(sq, "RGBA"), box


_font_warned: set[str] = set()


def font(name: str, size: int):
    """DejaVu from FONT_DIR, else the same face from the system font path, else Pillow's built-in font."""
    for cand in (FONT_DIR / name, Path(name)):
        try:
            return ImageFont.truetype(str(cand), size)
        except OSError:
            continue
    if name not in _font_warned:
        print(f"  WARNING: font {name} not found under {FONT_DIR} or on the font path; using Pillow's built-in font")
        _font_warned.add(name)
    try:
        return ImageFont.load_default(size=size)
    except TypeError:  # Pillow < 10.1 has no sized default font
        return ImageFont.load_default()


def wrap(text: str, f, width: float) -> list[str]:
    lines, cur = [], ""
    for word in text.split():
        cand = f"{cur} {word}".strip()
        if cur and f.getlength(cand) > width:
            lines.append(cur)
            cur = word
        else:
            cur = cand
    if cur:
        lines.append(cur)
    return lines


def fit(text: str, name: str, size: int, width: float):
    """The largest font size at or below `size` at which `text` fits in `width`."""
    f = font(name, size)
    while size > 12 and f.getlength(text) > width:
        size -= 1
        f = font(name, size)
    return f, size


def make_og(disc: Image.Image, question: str) -> Image.Image:
    """1200 x 630 card: the badge at the left, the question, an orange rule, the team line and the site URL at the right."""
    W, H = 1200, 630
    card = Image.new("RGBA", (W, H), CREAM + (255,))
    d = 470
    card.alpha_composite(disc.resize((d, d), LANCZOS), (70, (H - d) // 2))
    draw = ImageDraw.Draw(card)
    x0, col_w = 600, 530
    line_w = W - x0 - 40  # the single-line captions may run a little past the wrapped column
    size = 52
    while True:
        f_q = font("DejaVuSans-Bold.ttf", size)
        lines = wrap(question, f_q, col_w)
        if len(lines) <= 3 or size <= 30:
            break
        size -= 2
    lh = round(size * 1.18)
    f_team, team_size = fit(TEAM_LINE, "DejaVuSans.ttf", 30, line_w)
    f_url, url_size = fit(SITE_URL, "DejaVuSans.ttf", 26, line_w)
    gap_q, rule_h, gap_r, team_h, gap_t, url_h = 26, 6, 26, 36, 8, 32
    y = (H - (len(lines) * lh + gap_q + rule_h + gap_r + team_h + gap_t + url_h)) // 2
    for line in lines:
        draw.text((x0, y), line, font=f_q, fill=NAVY, anchor="la")
        y += lh
    y += gap_q
    draw.rectangle([x0, y, x0 + 140, y + rule_h - 1], fill=ORANGE)
    y += rule_h + gap_r
    draw.text((x0, y), TEAM_LINE, font=f_team, fill=INK2, anchor="la")
    y += team_h + gap_t
    draw.text((x0, y), SITE_URL, font=f_url, fill=INK2, anchor="la")
    print(f"  social card: question set at {size} px in {len(lines)} line(s), column {col_w} px; team line {team_size} px, "
          f"URL {url_size} px (fit to {line_w} px); badge {d} px")
    return card.convert("RGB")


def save(im: Image.Image, path: Path, **kw) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    im.save(path, **kw)
    n = path.stat().st_size
    warn = "   WARNING: above 400 KB" if n > SIZE_WARN else ""
    print(f"  {str(path.relative_to(ROOT)):36} {im.size[0]:>5} x {im.size[1]:<5} {n / 1000:8.1f} KB{warn}")


def main() -> int:
    if not SRC.exists():
        print(f"missing {SRC}")
        return 1
    with Image.open(SRC) as src:
        print(f"source {SRC.relative_to(ROOT)}: {src.size[0]} x {src.size[1]} {src.mode}")
        if src.size != (1600, 1800):
            print(f"  WARNING: expected 1600 x 1800; the expected badge geometry assumes that size")
        im = src.convert("RGBA")

    g = detect_badge(im)
    x0, y0, x1, y1 = g["bbox"]
    print(f"badge: centre ({g['cx']:.1f}, {g['cy']:.1f}), radius {g['r']:.1f} px; bbox x {x0}-{x1}, y {y0}-{y1} "
          f"({g['w']} x {g['h']}); {g['components']} components above the cream threshold, ring component {g['pixels']:,} px")
    offsets = {k: g[k] - EXPECTED[k] for k in ("cx", "cy", "r")}
    print(f"  expected centre ({EXPECTED['cx']:.0f}, {EXPECTED['cy']:.0f}) radius {EXPECTED['r']:.0f} ±{TOLERANCE}: "
          f"dx {offsets['cx']:+.1f}, dy {offsets['cy']:+.1f}, dr {offsets['r']:+.1f}")
    bad = [k for k, v in offsets.items() if abs(v) > TOLERANCE]
    if bad:
        print(f"  FAILED: {', '.join(bad)} outside the tolerance; check the source image")
        return 1

    disc, box = crop_disc(im, g)
    print(f"  crop {box} -> {disc.size[0]} px square ({MARGIN:.1%} margin); alpha mask r = {MASK_RADIUS} x "
          f"{disc.size[0] / 2:.1f} = {MASK_RADIUS * disc.size[0] / 2:.1f} px, feather {FEATHER} px")

    print("outputs:")
    for s in (16, 32, 64, 192, 512):
        save(disc.resize((s, s), LANCZOS), BRAND / f"badge-{s}.png", optimize=True)
    # iOS home-screen icon: opaque cream square, the disc inset so the rounded-corner mask does not clip the ring
    touch, inset = Image.new("RGBA", (180, 180), CREAM + (255,)), 164
    touch.alpha_composite(disc.resize((inset, inset), LANCZOS), ((180 - inset) // 2, (180 - inset) // 2))
    save(touch.convert("RGB"), BRAND / "badge-180.png", optimize=True)
    frames = [disc.resize((s, s), LANCZOS) for s in (48, 32, 16)]
    save(frames[0], FAVICON, format="ICO", sizes=[(48, 48), (32, 32), (16, 16)], append_images=frames[1:])
    with Image.open(FAVICON) as chk:
        print(f"    favicon frames: {sorted(chk.info.get('sizes', set()))}")
    save(im.convert("RGB").resize((800, 900), LANCZOS), BRAND / "ratpac-logo.png", optimize=True)

    if not HEADLINE.exists():
        print(f"missing {HEADLINE} (run `make site-data`): the social card needs its question")
        return 1
    question = json.loads(HEADLINE.read_text())["question"]
    save(make_og(disc, question), BRAND / "badge-og.png", optimize=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
