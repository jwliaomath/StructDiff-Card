from __future__ import annotations

import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageSequence

from structdiff_card.structure import read_structure

ROOT = Path(__file__).resolve().parents[1]
RESULT = json.loads((ROOT / "examples" / "result" / "result.json").read_text(encoding="utf-8"))
REFERENCE = read_structure(ROOT / "examples" / "result" / "reference.pdb")
MOBILE = read_structure(ROOT / "examples" / "result" / "mobile_aligned_to_reference.pdb")
OUTPUT = ROOT / "site" / "public" / "demo"

INK = "#07131b"
PAPER = "#f7faf9"
WHITE = "#edf7f5"
TEAL = "#42d5c5"
TEAL_DARK = "#0d9488"
CORAL = "#ff8a5c"
PURPLE = "#b77ce8"
YELLOW = "#ffcf56"
MUTED = "#8da3ac"
LINE = "#29434e"
RESIDUAL_STOPS = ["#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"]
RESIDUALS = {
    (row["mobile_chain"], row["mobile_residue"]): row["displacement"]
    for row in RESULT["residue_mapping"]
    if row["displacement"] is not None
}


def font(size: int, *, serif: bool = False, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = []
    if serif:
        candidates.extend([Path("C:/Windows/Fonts/georgiab.ttf" if bold else "C:/Windows/Fonts/georgia.ttf")])
    elif bold:
        candidates.extend([Path("C:/Windows/Fonts/arialbd.ttf")])
    else:
        candidates.extend([Path("C:/Windows/Fonts/arial.ttf")])
    candidates.extend(
        [
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf" if serif else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size)
    return ImageFont.load_default()


F10 = font(10, bold=True)
F12 = font(12)
F13 = font(13, bold=True)
F15 = font(15)
F18 = font(18, bold=True)
F24 = font(24, serif=True)
F34 = font(34, serif=True)
F46 = font(46, serif=True)
F62 = font(62, serif=True)


def fit_text(draw: ImageDraw.ImageDraw, text: str, box_width: int, base_size: int) -> ImageFont.FreeTypeFont:
    size = base_size
    while size > 12:
        chosen = font(size, serif=True)
        if draw.textbbox((0, 0), text, font=chosen)[2] <= box_width:
            return chosen
        size -= 2
    return font(size, serif=True)


def fade(value: float) -> float:
    return max(0.0, min(1.0, value))


def ease(value: float) -> float:
    value = fade(value)
    return value * value * (3 - 2 * value)


def residual_color(value: float) -> str:
    bounded = max(0.0, min(5.0, value))
    lower = min(4, int(bounded))
    fraction = bounded - lower
    channels = [
        tuple(int(RESIDUAL_STOPS[index][offset : offset + 2], 16) for offset in (1, 3, 5))
        for index in (lower, lower + 1)
    ]
    mixed = tuple(
        round(channels[0][index] + (channels[1][index] - channels[0][index]) * fraction)
        for index in range(3)
    )
    return "#" + "".join(f"{channel:02x}" for channel in mixed)


def draw_label(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, color: str = TEAL) -> None:
    draw.text(xy, text.upper(), fill=color, font=F10, spacing=2)


def rotate_project(
    xyz: tuple[float, float, float],
    center: tuple[float, float, float],
    angle: float,
    origin: tuple[float, float],
    scale: float,
) -> tuple[float, float, float]:
    x, y, z = (xyz[index] - center[index] for index in range(3))
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    x, z = cos_a * x + sin_a * z, -sin_a * x + cos_a * z
    tilt = -0.36
    cos_t, sin_t = math.cos(tilt), math.sin(tilt)
    y, z = cos_t * y - sin_t * z, sin_t * y + cos_t * z
    return origin[0] + scale * x, origin[1] - scale * y, z


def molecule_coordinates():
    points = []
    for residues in REFERENCE.chains.values():
        points.extend(residue.xyz for residue in residues)
    for residues in MOBILE.chains.values():
        points.extend(residue.xyz for residue in residues)
    center = tuple(sum(point[index] for point in points) / len(points) for index in range(3))
    spread = max(
        max(abs(point[index] - center[index]) for point in points)
        for index in range(3)
    )
    return center, 190 / spread


CENTER, MODEL_SCALE = molecule_coordinates()


def draw_molecule(draw: ImageDraw.ImageDraw, angle: float, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    origin = ((left + right) / 2, (top + bottom) / 2)
    records = []
    for chain, residues in REFERENCE.chains.items():
        projected = [rotate_project(residue.xyz, CENTER, angle, origin, MODEL_SCALE) for residue in residues]
        records.append((sum(item[2] for item in projected) / len(projected), "reference", chain, residues, projected))
    for chain, residues in MOBILE.chains.items():
        projected = [rotate_project(residue.xyz, CENTER, angle, origin, MODEL_SCALE) for residue in residues]
        kind = "mapped" if chain in {"C", "D"} else "unmatched"
        records.append((sum(item[2] for item in projected) / len(projected), kind, chain, residues, projected))
    for _, kind, chain, residues, projected in sorted(records, key=lambda item: item[0]):
        points = [(item[0], item[1]) for item in projected]
        if kind == "mapped":
            for index in range(1, len(points)):
                value = RESIDUALS.get((chain, residues[index].residue_id), 0.0)
                draw.line((points[index - 1], points[index]), fill=residual_color(value), width=4)
        else:
            color = "#9aabb2" if kind == "reference" else PURPLE
            draw.line(points, fill=color, width=2, joint="curve")

    threshold = RESULT["difference_summary"]["hotspot_threshold"]
    mobile_lookup = {
        (chain, residue.residue_id): residue
        for chain, residues in MOBILE.chains.items()
        for residue in residues
    }
    reference_lookup = {
        (chain, residue.residue_id): residue
        for chain, residues in REFERENCE.chains.items()
        for residue in residues
    }
    high_rows = sorted(
        (
            row
            for row in RESULT["residue_mapping"]
            if row["displacement"] is not None and row["displacement"] >= threshold
        ),
        key=lambda row: row["displacement"],
        reverse=True,
    )[:10]
    for row in high_rows:
        mobile_residue = mobile_lookup[(row["mobile_chain"], row["mobile_residue"])]
        reference_residue = reference_lookup[(row["reference_chain"], row["reference_residue"])]
        mobile_point = rotate_project(mobile_residue.xyz, CENTER, angle, origin, MODEL_SCALE)
        reference_point = rotate_project(reference_residue.xyz, CENTER, angle, origin, MODEL_SCALE)
        draw.line(
            ((mobile_point[0], mobile_point[1]), (reference_point[0], reference_point[1])),
            fill=YELLOW,
            width=2,
        )


def base_frame() -> tuple[Image.Image, ImageDraw.ImageDraw]:
    image = Image.new("RGB", (960, 540), INK)
    draw = ImageDraw.Draw(image)
    draw.rectangle((32, 26, 38, 61), fill=TEAL)
    draw.text((52, 28), "STRUCTDIFF CARD", fill=WHITE, font=F18)
    draw.text((52, 50), "MULTI-CHAIN STRUCTURE COMPARISON", fill=MUTED, font=F10)
    draw.line((32, 78, 928, 78), fill=LINE, width=1)
    return image, draw


def hero_frame(progress: float) -> Image.Image:
    image, _draw = base_frame()
    alpha = ease(progress)
    overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
    layer = ImageDraw.Draw(overlay)
    title_color = (237, 247, 245, int(255 * alpha))
    muted_color = (141, 163, 172, int(255 * alpha))
    layer.text((58, 148), "See where the", fill=title_color, font=F62)
    layer.text((58, 215), "assembly moves.", fill=(66, 213, 197, int(255 * alpha)), font=F62)
    layer.text((61, 314), "US-align → correspondence → contact, motion, and spatial differences", fill=muted_color, font=F15)
    layer.rectangle((60, 370, 285, 418), fill=(13, 148, 136, int(255 * alpha)))
    layer.text((82, 384), "COMPLETE EXAMPLE  ↓", fill=title_color, font=F13)
    image = Image.alpha_composite(image.convert("RGBA"), overlay).convert("RGB")
    return image


def overlay_frame(progress: float) -> Image.Image:
    image, draw = base_frame()
    draw_label(draw, (32, 96), "Interactive overlay")
    draw.text((32, 116), "2HHB → 1HHO", fill=WHITE, font=F34)
    draw.rectangle((30, 164, 615, 500), outline=LINE, width=1)
    draw_molecule(draw, angle=progress * math.tau * 0.9, box=(40, 168, 608, 493))
    draw.ellipse((49, 476, 57, 484), fill="#9aabb2")
    draw.text((64, 473), "reference ghost", fill=MUTED, font=F10)
    for offset in range(70):
        draw.line((160 + offset, 477, 160 + offset, 484), fill=residual_color(offset / 14), width=1)
    draw.text((238, 473), "0 → 5 Å", fill=MUTED, font=F10)
    draw.ellipse((306, 476, 314, 484), fill=YELLOW)
    draw.text((321, 473), "vectors", fill=MUTED, font=F10)
    draw.ellipse((380, 476, 388, 484), fill=PURPLE)
    draw.text((395, 473), "unmatched", fill=MUTED, font=F10)

    summary = RESULT["summary"]
    draw_label(draw, (650, 96), "Global result")
    metrics = [
        ("TM-SCORE · REF", f"{summary['tm_reference']:.3f}"),
        ("RMSD", f"{summary['rmsd']:.2f} Å"),
        ("ALIGNED", f"{summary['aligned_length']} / {summary['reference_length']}"),
    ]
    y = 128
    for label, value in metrics:
        draw.text((650, y), label, fill=MUTED, font=F10)
        draw.text((650, y + 18), value, fill=WHITE, font=F34)
        y += 76
    draw.line((650, 363, 922, 363), fill=LINE)
    draw_label(draw, (650, 381), "Chain mapping")
    for index, row in enumerate(RESULT["chain_alignments"]):
        y = 410 + index * 38
        draw.rectangle((650, y, 680, y + 27), fill="#0e5b57")
        draw.text((660, y + 5), row["mobile_chain"], fill=TEAL, font=F13)
        draw.text((693, y + 5), "→", fill=MUTED, font=F13)
        draw.rectangle((718, y, 748, y + 27), fill="#663423")
        draw.text((728, y + 5), row["reference_chain"], fill=CORAL, font=F13)
        draw.text((770, y + 6), f"TM {row['tm_reference']:.3f}", fill=WHITE, font=F12)
    return image


def plot_frame(progress: float) -> Image.Image:
    image, draw = base_frame()
    draw_label(draw, (32, 96), "Per-residue displacement")
    draw.text((32, 116), "One global transform. No hidden chain refit.", fill=WHITE, font=F24)
    plot = (72, 184, 924, 412)
    rows = [row for row in RESULT["residue_mapping"] if row["displacement"] is not None]
    max_value = 5.8
    for tick in range(6):
        y = plot[3] - (tick / max_value) * (plot[3] - plot[1])
        draw.line((plot[0], y, plot[2], y), fill=CORAL if tick == 5 else LINE, width=1)
        draw.text((45, y - 6), str(tick), fill=MUTED, font=F10)
    visible = max(2, int(len(rows) * ease(progress)))
    points = []
    for index, row in enumerate(rows[:visible]):
        x = plot[0] + (index / max(1, len(rows) - 1)) * (plot[2] - plot[0])
        y = plot[3] - (row["displacement"] / max_value) * (plot[3] - plot[1])
        points.append((x, y))
    if len(points) > 1:
        draw.line(points, fill=TEAL, width=3, joint="curve")
    draw.text((72, 443), "matched residue pairs in C→A, D→B mapping order", fill=MUTED, font=F12)
    draw.rectangle((72, 475, 279, 505), fill="#1b152b")
    draw.text((83, 483), "A, B  ·  UNMATCHED CHAINS", fill="#d8b7f0", font=F10)
    return image


def diagnostic_frame(progress: float) -> Image.Image:
    image, draw = base_frame()
    alpha = ease(progress)
    draw_label(draw, (32, 96), "Structural explanation")
    draw.text((32, 116), "What changed, and how?", fill=WHITE, font=F34)
    summary = RESULT["contact_summary"]
    cards = [(32, 185, 309, 485), (341, 185, 618, 485), (650, 185, 927, 485)]
    for card in cards:
        draw.rectangle(card, fill="#0b1b23", outline=LINE)

    draw.text((52, 207), "01  CONTACTS", fill=TEAL, font=F10)
    draw.text((52, 242), f"{summary['gained_contacts']}", fill=TEAL, font=F46)
    draw.text((126, 261), "gained", fill=MUTED, font=F13)
    draw.text((52, 315), f"{summary['lost_contacts']}", fill=CORAL, font=F46)
    draw.text((126, 334), "lost", fill=MUTED, font=F13)
    draw.text((52, 412), "8 Å contact · 1 Å margin", fill=MUTED, font=F12)

    draw.text((361, 207), "02  CHAIN MOTION", fill=TEAL, font=F10)
    for index, row in enumerate(RESULT["motion_decomposition"]):
        y = 252 + index * 96
        draw.text((361, y), f"{row['mobile_chain']} → {row['reference_chain']}", fill=WHITE, font=F13)
        rigid_width = int(175 * row["rigid_body_rmsd"])
        internal_width = int(175 * row["internal_rmsd"])
        draw.rectangle((361, y + 27, 361 + rigid_width, y + 38), fill="#287fb8")
        draw.rectangle((361, y + 44, 361 + internal_width, y + 55), fill=CORAL)
        draw.text((545, y + 23), f"{row['rigid_body_rmsd']:.2f}", fill="#73b7e3", font=F10)
        draw.text((545, y + 40), f"{row['internal_rmsd']:.2f}", fill="#ff9a7b", font=F10)

    draw.text((670, 207), "03  3D PATCHES", fill=TEAL, font=F10)
    draw.text((670, 244), str(len(RESULT["spatial_patches"])), fill=YELLOW, font=F46)
    draw.text((735, 263), "spatial clusters", fill=MUTED, font=F13)
    for index, patch in enumerate(RESULT["spatial_patches"][:3]):
        y = 330 + index * 48
        draw.ellipse((670, y, 694, y + 24), fill=YELLOW)
        draw.text((677, y + 5), f"P{patch['rank']}", fill=INK, font=F10)
        draw.text(
            (706, y + 2),
            f"{patch['member_count']} residues · peak {patch['peak_displacement']:.2f} Å",
            fill=WHITE,
            font=F12,
        )
    if alpha < 1:
        veil = Image.new("RGBA", image.size, (7, 19, 27, int(255 * (1 - alpha))))
        image = Image.alpha_composite(image.convert("RGBA"), veil).convert("RGB")
    return image


def export_frame(progress: float) -> Image.Image:
    image, draw = base_frame()
    draw_label(draw, (32, 96), "Publication-ready output")
    draw.text((32, 116), "Download the evidence, not just a screenshot.", fill=WHITE, font=F34)
    card = (32, 182, 594, 492)
    draw.rectangle(card, fill="#0b1b23", outline=LINE)
    draw.rectangle((50, 202, 55, 245), fill=TEAL)
    draw.text((70, 204), "2HHB.PDB → 1HHO.PDB", fill=WHITE, font=F24)
    draw.text((70, 234), "MULTI-CHAIN ALIGNMENT", fill=MUTED, font=F10)
    draw.line((50, 270, 576, 270), fill=LINE)
    draw.text((52, 292), "TM-REF", fill=MUTED, font=F10)
    draw.text((52, 310), "0.982", fill=TEAL, font=F34)
    draw.text((235, 292), "RMSD", fill=MUTED, font=F10)
    draw.text((235, 310), "0.89 Å", fill=WHITE, font=F34)
    draw.text((423, 292), "ALIGNED", fill=MUTED, font=F10)
    draw.text((423, 310), "287", fill=WHITE, font=F34)
    draw.text((52, 380), "CHAIN MAP", fill=MUTED, font=F10)
    draw.text((52, 398), "C > A  ·  D > B", fill=WHITE, font=F18)
    draw.text((52, 442), "UNMATCHED", fill=MUTED, font=F10)
    draw.text((148, 440), "A, B", fill=PURPLE, font=F15)

    right_x = 635
    draw_label(draw, (right_x, 188), "Run locally")
    commands = [
        ("CSV", "residue + contact mapping"),
        ("CSV", "motion + spatial patches"),
        ("SVG", "plot + card"),
        ("PDB", "aligned + difference coded"),
        ("JSON", "full provenance"),
    ]
    for index, (kind, label) in enumerate(commands):
        y = 210 + index * 37
        draw.rectangle((right_x, y, right_x + 48, y + 27), fill="#15303a")
        draw.text((right_x + 9, y + 6), kind, fill=TEAL, font=F10)
        draw.text((right_x + 62, y + 6), label, fill=WHITE, font=F12)
    reveal = ease(progress)
    draw.rectangle((right_x, 408, 921, 477), fill="#020a0e", outline=LINE)
    command = "docker run --rm -p 8501:8501\nghcr.io/jwliaomath/structdiff-card:latest"
    chars = int(len(command) * reveal)
    draw.multiline_text((right_x + 14, 421), command[:chars], fill="#c8f6ef", font=F10, spacing=7)
    return image


def make_gif() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for index in range(100):
        if index < 14:
            frame = hero_frame(index / 11)
        elif index < 43:
            frame = overlay_frame((index - 14) / 29)
        elif index < 63:
            frame = plot_frame((index - 43) / 20)
        elif index < 81:
            frame = diagnostic_frame((index - 63) / 12)
        else:
            frame = export_frame((index - 81) / 15)
        frames.append(frame.quantize(colors=192, method=Image.Quantize.MEDIANCUT))
    frames[0].save(
        OUTPUT / "structdiff-demo.gif",
        save_all=True,
        append_images=frames[1:],
        duration=150,
        loop=0,
        optimize=True,
        disposal=2,
    )


def make_og_card() -> None:
    image = Image.new("RGB", (1200, 630), INK)
    draw = ImageDraw.Draw(image)
    draw.rectangle((62, 64, 70, 134), fill=TEAL)
    draw.text((94, 68), "STRUCTDIFF CARD", fill=WHITE, font=font(30, bold=True))
    draw.text((95, 112), "LOCAL-FIRST MULTI-CHAIN STRUCTURE COMPARISON", fill=MUTED, font=font(13, bold=True))
    draw.line((62, 166, 1138, 166), fill=LINE, width=2)
    title = "See where the assembly moves."
    draw.text((62, 210), title, fill=WHITE, font=fit_text(draw, title, 770, 64))
    draw.text(
        (65, 302),
        "Contact changes · chain motion · spatial patches · publication exports",
        fill=MUTED,
        font=font(19),
    )
    draw.rectangle((62, 380, 290, 510), outline=LINE, width=2)
    draw.text((84, 402), "TM-SCORE · REF", fill=MUTED, font=font(13, bold=True))
    draw.text((84, 438), "0.982", fill=TEAL, font=font(48, serif=True))
    draw.rectangle((314, 380, 542, 510), outline=LINE, width=2)
    draw.text((336, 402), "GLOBAL RMSD", fill=MUTED, font=font(13, bold=True))
    draw.text((336, 438), "0.89 Å", fill=WHITE, font=font(48, serif=True))
    draw.rectangle((566, 380, 794, 510), outline=LINE, width=2)
    draw.text((588, 402), "CHAIN MAP", fill=MUTED, font=font(13, bold=True))
    draw.text((588, 448), "C > A · D > B", fill=WHITE, font=font(25, bold=True))
    draw_molecule(draw, angle=0.65, box=(810, 180, 1170, 565))
    draw.text((1138, 588), "github.com/jwliaomath/StructDiff-Card", fill=MUTED, font=font(13), anchor="ra")
    image.save(OUTPUT / "og-card.png", optimize=True)


if __name__ == "__main__":
    make_gif()
    make_og_card()
    gif = Image.open(OUTPUT / "structdiff-demo.gif")
    duration = sum(frame.info.get("duration", 0) for frame in ImageSequence.Iterator(gif))
    print(f"GIF frames={gif.n_frames}, duration={duration} ms")
