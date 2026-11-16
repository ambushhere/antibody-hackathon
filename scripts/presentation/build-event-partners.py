"""Compose official event partner artwork without redrawing the logos."""
import base64
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "docs/assets/brands"
OUT = ROOT / "docs/assets/event-partners.svg"

PARTNERS = [
    ("apheris.png", "Apheris", "Main sponsor"),
    ("openfold.png", "OpenFold", "Silver sponsor"),
    ("nebius.png", "Nebius", "GPU resources"),
    ("tenstorrent-aiand.png", "Tenstorrent × ai&", "JapanFold credits"),
    ("amass.svg", "Amass", "Scientific intelligence"),
    ("adaptyv.png", "Adaptyv", "Laboratory validation"),
    ("anthropic.png", "Anthropic", "AI credits"),
    ("bi-startup-village.png", "BI Startup Village", "Venue"),
    ("oslo-cancer-cluster.png", "Oslo Cancer Cluster", "Event support"),
    ("girls-who-ml.png", "Girls Who ML", "Community and event design"),
    ("tekna-biotech.png", "Tekna Biotek", "Professional network"),
    ("yogutide.png", "Yogutide", "Hospitality"),
]


def artwork(name, x, y, width, height):
    mime = "image/svg+xml" if name.endswith(".svg") else "image/png"
    raw = base64.b64encode((ASSETS / name).read_bytes()).decode()
    return (
        f'<image x="{x}" y="{y}" width="{width}" height="{height}" '
        f'preserveAspectRatio="xMidYMid meet" href="data:{mime};base64,{raw}"/>'
    )


def build():
    parts = ['''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="590" viewBox="0 0 1200 590" role="img" aria-labelledby="title desc">
<title id="title">AI × Antibodies, Oslo 2026: organizer and event partners</title>
<desc id="desc">Organized by Nucleate Norway. Official logos acknowledge event partners and their roles.</desc>
<rect width="1200" height="590" rx="20" fill="#101d2c"/>
<path d="M40 146H1160" stroke="#334457"/>
<g font-family="Arial,Helvetica,sans-serif">
<text x="40" y="42" font-size="14" letter-spacing="2" fill="#8fb2c4">BUILT DURING THE HACKATHON</text>
<text x="38" y="93" font-size="45" font-weight="700" fill="#f6f3e9">AI × Antibodies</text>
<text x="40" y="123" font-size="18" fill="#b7c9d7">Oslo · 14 November 2026</text>
<text x="992" y="36" font-size="12" letter-spacing="1.3" fill="#b7c9d7">ORGANIZED BY</text>
''', artwork("nucleate-norway.svg", 980, 48, 170, 82)]
    for i, (name, brand, role) in enumerate(PARTNERS):
        x = 40 + (i % 4) * 284
        y = 166 + (i // 4) * 118
        parts += [
            f'<g><title>{escape(brand)} — {escape(role)}</title>',
            f'<rect x="{x}" y="{y}" width="268" height="102" rx="10" fill="#162537" stroke="#273c50"/>',
        ]
        if name == "girls-who-ml.png":
            parts.append(f'<rect x="{x+20}" y="{y+10}" width="228" height="55" rx="4" fill="#f4f1e8"/>')
        parts += [artwork(name, x+24, y+12, 220, 50),
                  f'<text x="{x+134}" y="{y+84}" text-anchor="middle" font-size="13" fill="#c3d1dd">{escape(role)}</text></g>']
    parts += ['<text x="40" y="558" font-size="13" fill="#a1b4c3">Event partners · original artwork from the official event website</text></g></svg>']
    OUT.write_text("\n".join(parts))
    print(OUT)


if __name__ == "__main__":
    build()
