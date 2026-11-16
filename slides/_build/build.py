"""Build self-contained slide decks.

Each  <name>.src.html  in this folder is turned into  ../<dayN>/<name>.html
with d3, the GBSG2 data and the GIMM logos inlined, so the result works offline.
The day folder is taken from the file name prefix (Day1_..., Day2_...).

Usage:  python slides/_build/build.py
"""
import base64, json, pathlib, re

HERE = pathlib.Path(__file__).resolve().parent
d3 = (HERE / "d3.min.js").read_text(encoding="utf-8")
data = json.dumps(json.loads((HERE / "gbsg2.json").read_text()), separators=(",", ":"))
logo = lambda f: "data:image/png;base64," + base64.b64encode((HERE / f).read_bytes()).decode()

for src in sorted(HERE.glob("*.src.html")):
    name = src.name.replace(".src.html", "")
    day = re.match(r"(Day\d+)_", name).group(1).lower()
    html = (src.read_text(encoding="utf-8")
            .replace("/*D3*/", d3).replace("__DATA__", data)
            .replace("__LOGO_NEG__", logo("logo_neg.png"))
            .replace("__LOGO_POS__", logo("logo_pos.png")))
    out = HERE.parent / day / f"{name}.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"built {out.relative_to(HERE.parent.parent)}  ({len(html)//1024} KB)")
