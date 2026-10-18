"""
Rich-text documentation viewer for the XPS Thickness Suite.

The documentation lives once, as Markdown with LaTeX maths in docs/*.md (GitHub renders it natively). This module
shows the same files inside the app: headings, lists, tables, images and links in a Tk Text widget, with the
equations typeset by matplotlib's mathtext (readable plain-text equations if matplotlib is missing).

It deliberately supports only a SUBSET of Markdown (see parse_markdown). parse_markdown(..., strict=True) raises on
anything outside that subset, and tests/test_docs.py runs every page through it, so the docs cannot drift from what
the app can display. Supported:
  blocks   # / ## / ### headings; paragraphs; "- " / "* " / "1. " lists (two levels); "> " quotes; ``` code fences;
           ```math fences (display equations); pipe tables (plain text cells); ![alt](path) images on their own
           line; --- rules (after a blank line); <!-- comments -->
  inline   **bold**, *italic*, `code`, [text](url), $math$, and backslash escapes. Bare http(s) URLs are links.
  maths    only what BOTH GitHub (KaTeX) and matplotlib mathtext understand: \\frac, \\sqrt, \\left \\right, \\ln \\exp
           \\sinh \\cos ..., Greek letters, sub/superscripts, \\cdot \\times \\approx \\le \\ge \\pm \\infty, \\mathrm{}.

The parsing half is pure Python (no Tk) so it can be tested without a display.
"""
import base64
import io
import os
import re
import webbrowser
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

PAGES = [("User guide", "user-guide.md"), ("Methods", "methods.md"), ("Presets & data", "presets.md"),
         ("References", "references.md")]
GITHUB_DOCS_URL = "https://github.com/aktiv8/XPS-Thickness-Suite/tree/main/docs"


# ---------------------------------------------------------------------------
# Markdown parsing (pure Python)
# ---------------------------------------------------------------------------

@dataclass
class Block:
    kind: str                       # heading para list quote code math table image rule
    text: str = ""
    level: int = 0                  # heading level
    lang: str = ""
    items: List[Tuple[int, str, str]] = field(default_factory=list)      # list: (level, marker, text)
    header: List[str] = field(default_factory=list)
    rows: List[List[str]] = field(default_factory=list)
    alt: str = ""
    path: str = ""


@dataclass(frozen=True)
class Span:
    kind: str                       # text code math link
    text: str
    bold: bool = False
    italic: bool = False
    url: str = ""


_FENCE = re.compile(r"^```(\w*)\s*$")
_HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*#*\s*$")
_LIST = re.compile(r"^(\s*)([-*+]|\d+[.)])\s+(.*)$")
_IMAGE = re.compile(r"^!\[([^\]]*)\]\(([^)\s]+)\)\s*$")
_RULE = re.compile(r"^\s*(-{3,}|\*{3,}|_{3,})\s*$")
_TABLE_SEP = re.compile(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)*\|?\s*$")
_LINK = re.compile(r"\[([^\]]+)\]\(([^)\s]+)\)")
_URL = re.compile(r"https?://[^\s<>)\]]+")
_HTML = re.compile(r"<[A-Za-z/!]")
_UNDERSCORE_EM = re.compile(r"(?<![A-Za-z0-9_])_[^_\s](?:[^_]*[^_\s])?_(?![A-Za-z0-9_])")


def _is_block_start(lines, i):
    line = lines[i]
    return bool(_FENCE.match(line) or _HEADING.match(line) or _LIST.match(line) or _IMAGE.match(line)
                or line.startswith(">") or _RULE.match(line) or _is_table_start(lines, i)
                or line.lstrip().startswith("<!--"))


def _is_table_start(lines, i):
    return "|" in lines[i] and i + 1 < len(lines) and bool(_TABLE_SEP.match(lines[i + 1])) and "-" in lines[i + 1]


def _split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


def parse_markdown(text: str, strict: bool = False) -> List[Block]:
    """Parse the supported Markdown subset into blocks. strict=True raises ValueError on anything else."""
    lines = text.replace("\r\n", "\n").split("\n")
    blocks: List[Block] = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        m = _FENCE.match(line)
        if m:
            body, j = [], i + 1
            while j < n and not lines[j].startswith("```"):
                body.append(lines[j])
                j += 1
            if j >= n:
                raise ValueError(f"Unclosed code fence starting at line {i + 1}")
            blocks.append(Block("math", text="\n".join(body).strip()) if m.group(1) == "math"
                          else Block("code", lang=m.group(1), text="\n".join(body)))
            i = j + 1
            continue
        if line.lstrip().startswith("<!--"):
            while i < n and "-->" not in lines[i]:
                i += 1
            i += 1
            continue
        m = _HEADING.match(line)
        if m:
            level = len(m.group(1))
            if level > 3 and strict:
                raise ValueError(f"Headings deeper than ### are not supported (line {i + 1}: {line!r})")
            blocks.append(Block("heading", text=m.group(2), level=min(level, 3)))
            i += 1
            continue
        if _RULE.match(line):
            if strict and i > 0 and lines[i - 1].strip():
                raise ValueError(f"A horizontal rule needs a blank line before it (line {i + 1})")
            blocks.append(Block("rule"))
            i += 1
            continue
        m = _IMAGE.match(line)
        if m:
            blocks.append(Block("image", alt=m.group(1), path=m.group(2)))
            i += 1
            continue
        if line.startswith(">"):
            quote = []
            while i < n and lines[i].startswith(">"):
                quote.append(re.sub(r"^>\s?", "", lines[i]))
                i += 1
            blocks.append(Block("quote", text=" ".join(q.strip() for q in quote if q.strip())))
            continue
        if _is_table_start(lines, i):
            header = _split_row(line)
            rows, i = [], i + 2
            while i < n and "|" in lines[i] and lines[i].strip():
                row = _split_row(lines[i])
                if strict and len(row) != len(header):
                    raise ValueError(f"Table row has {len(row)} cells, header has {len(header)}: {lines[i]!r}")
                rows.append(row)
                i += 1
            if strict and any("$" in c for c in header + [c for r in rows for c in r]):
                raise ValueError("Maths is not supported inside table cells (use plain text or Unicode there).")
            blocks.append(Block("table", header=header, rows=rows))
            continue
        m = _LIST.match(line)
        if m:
            items: List[Tuple[int, str, str]] = []
            ordered = m.group(2)[0].isdigit()          # a switch between bullets and numbers starts a new list
            while i < n:
                m = _LIST.match(lines[i])
                if m:
                    level = len(m.group(1).expandtabs(2)) // 2
                    if level > 1 and strict:
                        raise ValueError(f"Lists nested deeper than two levels are not supported (line {i + 1})")
                    marker = m.group(2)
                    if level == 0 and marker[0].isdigit() != ordered:
                        break
                    items.append((min(level, 1), "•" if marker in "-*+" else marker, m.group(3)))
                    i += 1
                elif lines[i].strip() and lines[i][0] in " \t" and items:       # continuation of the previous item
                    lvl, mk, tx = items[-1]
                    items[-1] = (lvl, mk, tx + " " + lines[i].strip())
                    i += 1
                elif not lines[i].strip() and i + 1 < n and _LIST.match(lines[i + 1]):
                    i += 1                                                      # a blank line inside a loose list
                else:
                    break
            blocks.append(Block("list", items=items))
            continue
        para = [line.strip()]
        i += 1
        while i < n and lines[i].strip() and not _is_block_start(lines, i):
            para.append(lines[i].strip())
            i += 1
        blocks.append(Block("para", text=" ".join(para)))
    return blocks


def parse_inline(text: str, strict: bool = False) -> List[Span]:
    """Split inline Markdown into spans (text, code, math, link) carrying bold/italic flags."""
    if strict:
        if _HTML.search(text):
            raise ValueError(f"Raw HTML is not supported: {text!r}")
        if "~~" in text:
            raise ValueError(f"Strikethrough is not supported: {text!r}")
        stripped = re.sub(r"`[^`]*`|\$[^$]*\$", "", text)
        if _UNDERSCORE_EM.search(stripped):
            raise ValueError(f"_emphasis_ is not supported (use *italic*): {text!r}")
    return _inline(text, False, False, strict)


def _autolink(text, bold, italic):
    spans, pos = [], 0
    for m in _URL.finditer(text):
        url = m.group(0).rstrip(".,;:")
        if m.start() > pos:
            spans.append(Span("text", text[pos:m.start()], bold, italic))
        spans.append(Span("link", url, bold, italic, url=url))
        pos = m.start() + len(url)
    if pos < len(text):
        spans.append(Span("text", text[pos:], bold, italic))
    return spans


def _inline(text, bold, italic, strict):
    spans: List[Span] = []
    buf: List[str] = []
    i, n = 0, len(text)

    def flush():
        if buf:
            spans.extend(_autolink("".join(buf), bold, italic))
            buf.clear()

    def bad(msg):
        if strict:
            raise ValueError(f"{msg} in: {text!r}")

    while i < n:
        c = text[i]
        if c == "\\" and i + 1 < n and text[i + 1] in "\\`*_{}[]()#+-.!$|<>~":
            buf.append(text[i + 1])
            i += 2
            continue
        if c == "`":
            j = text.find("`", i + 1)
            if j < 0:
                bad("Unmatched `")
                buf.append(c)
                i += 1
                continue
            flush()
            spans.append(Span("code", text[i + 1:j], bold, italic))
            i = j + 1
            continue
        if c == "$":
            j = i + 1
            while j < n and not (text[j] == "$" and text[j - 1] != "\\"):
                j += 1
            if j >= n or j == i + 1 or text[i + 1].isspace() or text[j - 1].isspace():
                bad("Unmatched or malformed $...$")
                buf.append(c)
                i += 1
                continue
            flush()
            spans.append(Span("math", text[i + 1:j], bold, italic))
            i = j + 1
            continue
        if c == "[":
            m = _LINK.match(text, i)
            if m:
                flush()
                spans.append(Span("link", m.group(1), bold, italic, url=m.group(2)))
                i = m.end()
                continue
        if text.startswith("**", i):
            j = text.find("**", i + 2)
            if j < 0 or text[i + 2:i + 3].isspace():
                bad("Unmatched **")
                buf.append("**")
                i += 2
                continue
            flush()
            spans.extend(_inline(text[i + 2:j], True, italic, strict))
            i = j + 2
            continue
        if c == "*":
            j = i + 1
            while j < n and not (text[j] == "*" and text[j - 1] != "\\" and not text.startswith("**", j)):
                j += 1
            if j >= n or text[i + 1:i + 2].isspace() or text[j - 1].isspace():
                bad("Unmatched or malformed *italic*")
                buf.append(c)
                i += 1
                continue
            flush()
            spans.extend(_inline(text[i + 1:j], bold, True, strict))
            i = j + 1
            continue
        buf.append(c)
        i += 1
    flush()
    return spans


def slugify(heading: str) -> str:
    """GitHub-style anchor for a heading."""
    plain = re.sub(r"[`*]|\$[^$]*\$", "", heading)
    return re.sub(r"[^\w\- ]", "", plain.lower()).replace(" ", "-")


def heading_slugs(blocks: List[Block]) -> List[str]:
    """Anchors of all headings, with GitHub's -1, -2 suffixes for duplicates."""
    seen, out = {}, []
    for b in blocks:
        if b.kind == "heading":
            s = slugify(b.text)
            k = seen.get(s, 0)
            seen[s] = k + 1
            out.append(s if k == 0 else f"{s}-{k}")
    return out


def iter_inline_texts(blocks: List[Block]):
    """Every inline-markup string in the blocks (not code/math blocks)."""
    for b in blocks:
        if b.kind in ("heading", "para", "quote"):
            yield b.text
        elif b.kind == "list":
            for _, _, t in b.items:
                yield t
        elif b.kind == "table":
            yield from b.header
            for r in b.rows:
                yield from r


def normalise_math(expr: str) -> str:
    """One line, single spaces: GitHub accepts multi-line maths, matplotlib's mathtext does not."""
    return re.sub(r"\s+", " ", expr).strip()


def iter_math(blocks: List[Block]):
    """Every LaTeX expression in the blocks (```math blocks and inline $...$), whitespace-normalised."""
    for b in blocks:
        if b.kind == "math":
            yield normalise_math(b.text)
    for t in iter_inline_texts(blocks):
        for s in parse_inline(t):
            if s.kind == "math":
                yield normalise_math(s.text)


# ---------------------------------------------------------------------------
# Plain-text fallback for equations (used when matplotlib is not installed)
# ---------------------------------------------------------------------------

_GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "eta": "η", "theta": "θ",
          "lambda": "λ", "mu": "μ", "nu": "ν", "pi": "π", "rho": "ρ", "sigma": "σ", "tau": "τ", "phi": "φ",
          "chi": "χ", "psi": "ψ", "omega": "ω", "xi": "ξ", "Delta": "Δ", "Sigma": "Σ", "Omega": "Ω",
          "cdot": "·", "times": "×", "approx": "≈", "le": "≤", "ge": "≥", "pm": "±", "infty": "∞", "sim": "~",
          "rightarrow": "→", "to": "→", "sum": "Σ", "int": "∫", "partial": "∂", "leq": "≤", "geq": "≥",
          "lesssim": "≲", "gtrsim": "≳", "langle": "⟨", "rangle": "⟩", "propto": "∝", "ldots": "…",
          "varepsilon": "ε", "kappa": "κ", "zeta": "ζ", "Lambda": "Λ", "Gamma": "Γ", "Theta": "Θ", "Pi": "Π"}
_SUP = str.maketrans("0123456789+-=()n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿ")
_SUB = str.maketrans("0123456789+-=()", "₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎")


def _group(s, i):
    """Return (content, next index) of a {...} group or a single character at s[i:]."""
    if i < len(s) and s[i] == "{":
        depth, j = 0, i
        while j < len(s):
            depth += (s[j] == "{") - (s[j] == "}")
            if depth == 0:
                return s[i + 1:j], j + 1
            j += 1
        return s[i + 1:], len(s)
    return (s[i], i + 1) if i < len(s) else ("", i)


def latex_to_text(expr: str) -> str:
    """A readable Unicode approximation of a LaTeX expression (the no-matplotlib fallback)."""
    out, i, s = [], 0, expr
    while i < len(s):
        c = s[i]
        if c == "\\":
            m = re.match(r"\\([A-Za-z]+|.)", s[i:])
            name = m.group(1)
            i += len(m.group(0))
            if name == "frac":
                a, i = _group(s, i)
                b, i = _group(s, i)
                out.append(f"({latex_to_text(a)})/({latex_to_text(b)})")
            elif name == "sqrt":
                a, i = _group(s, i)
                out.append(f"√({latex_to_text(a)})")
            elif name in ("mathrm", "text", "mathbf", "mathit"):
                a, i = _group(s, i)
                out.append(latex_to_text(a))
            elif name in ("left", "right"):
                pass
            elif name in _GREEK:
                out.append(_GREEK[name])
            elif name in (",", ";", ":", "!", " ", "quad", "qquad"):
                out.append(" " if name not in ("!",) else "")
            elif name in ("ln", "exp", "sinh", "cosh", "sin", "cos", "tan", "log", "max", "min"):
                out.append(name)
            else:
                out.append(name)
        elif c in "^_":
            a, i2 = _group(s, i + 1)
            i = i2
            a = latex_to_text(a)
            table = _SUP if c == "^" else _SUB
            conv = a.translate(table)
            out.append(conv if all(ch in "0123456789+-=()n" for ch in a) and c == "^" or
                       all(ch in "0123456789+-=()" for ch in a) and c == "_" else (f"^({a})" if c == "^" else f"_{a}"))
        elif c in "{}":
            i += 1
        else:
            out.append(c)
            i += 1
    return re.sub(r"\s+", " ", "".join(out)).strip()


# ---------------------------------------------------------------------------
# Tk viewer
# ---------------------------------------------------------------------------

def docs_dir() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "docs")


class MathRenderer:
    """Typeset LaTeX with matplotlib mathtext into cached Tk PhotoImages; returns None if it cannot."""

    def __init__(self):
        self.cache = {}
        try:
            from matplotlib import mathtext
            from matplotlib.font_manager import FontProperties
            import matplotlib
            self._mathtext, self._fp, self._mpl = mathtext, FontProperties, matplotlib
            self.available = True
        except Exception:
            self.available = False

    def image(self, expr: str, display: bool, color: str = "#1f2937"):
        import tkinter as tk
        expr = normalise_math(expr)
        key = (expr, display, color)
        if key in self.cache:
            return self.cache[key]
        img = None
        if self.available:
            try:
                buf = io.BytesIO()
                with self._mpl.rc_context({"mathtext.fontset": "stix"}):
                    self._mathtext.math_to_image(f"${expr}$", buf, prop=self._fp(size=15 if display else 11.5),
                                                 dpi=220, format="png", color=color)
                try:                     # render at 2x and shrink with a smooth filter: crisper small exponents
                    from PIL import Image, ImageTk
                    big = Image.open(io.BytesIO(buf.getvalue())).convert("RGBA")
                    img = ImageTk.PhotoImage(big.resize((max(1, big.width // 2), max(1, big.height // 2)),
                                                        Image.LANCZOS))
                except ImportError:       # no Pillow: Tk's own 2x subsample
                    img = tk.PhotoImage(data=base64.b64encode(buf.getvalue())).subsample(2)
            except Exception:
                img = None
        self.cache[key] = img
        return img


def build_viewer(parent, strict: bool = False, math: Optional[MathRenderer] = None, start: Optional[str] = None):
    """Create the documentation viewer frame (page list + rich text). Imported lazily so the module stays Tk-free."""
    import tkinter as tk
    import tkinter.font as tkfont
    from tkinter import ttk
    from gui_common import ACCENT, CARD, INK, MUTED

    class DocView(ttk.Frame):
        BODY, MONO = "Segoe UI", "Consolas"

        def __init__(self):
            super().__init__(parent)
            self.strict = strict
            self.math = math or MathRenderer()
            self.dir = docs_dir()
            self.page = None
            self._keep = []             # PhotoImages / fonts must stay referenced
            self._anchors = {}
            self.columnconfigure(1, weight=1)
            self.rowconfigure(0, weight=1)

            side = ttk.Frame(self)
            side.grid(row=0, column=0, sticky="ns", padx=(0, 8), pady=4)
            ttk.Label(side, text="Documentation", style="RowHead.TLabel").pack(anchor="w", pady=(2, 4))
            self.listbox = tk.Listbox(side, exportselection=False, activestyle="none", width=17, height=len(PAGES),
                                      font=(self.BODY, 10, "bold"), relief="flat", highlightthickness=1,
                                      highlightbackground="#d1d5db", selectbackground=ACCENT, selectforeground="white",
                                      bg=CARD, fg=INK)
            for title, _ in PAGES:
                self.listbox.insert("end", " " + title)
            self.listbox.pack(fill="x")
            self.listbox.bind("<<ListboxSelect>>", lambda e: self._on_select())
            ttk.Button(side, text="Open on GitHub", command=lambda: webbrowser.open(GITHUB_DOCS_URL)).pack(
                fill="x", pady=(10, 0))
            note = ("Equations typeset" if self.math.available else "Install matplotlib to typeset equations")
            ttk.Label(side, text=note, style="Muted.TLabel", wraplength=130, justify="left").pack(anchor="w", pady=8)

            body = ttk.Frame(self)
            body.grid(row=0, column=1, sticky="nsew", pady=4)
            body.columnconfigure(0, weight=1)
            body.rowconfigure(0, weight=1)
            self.text = tk.Text(body, wrap="word", relief="flat", bg=CARD, fg=INK, padx=22, pady=14,
                                highlightthickness=1, highlightbackground="#d1d5db", cursor="arrow",
                                font=(self.BODY, 10), spacing1=2, spacing3=4, state="disabled")
            scroll = ttk.Scrollbar(body, orient="vertical", command=self.text.yview)
            self.text.configure(yscrollcommand=scroll.set)
            self.text.grid(row=0, column=0, sticky="nsew")
            scroll.grid(row=0, column=1, sticky="ns")
            self._configure_tags()
            self.show(start or PAGES[0][1])

        # --- tags --------------------------------------------------------------------------------
        def _font(self, **kw):
            f = tkfont.Font(family=kw.pop("family", self.BODY), **kw)
            self._keep.append(f)
            return f

        def _configure_tags(self):
            t = self.text
            t.tag_configure("h1", font=self._font(size=19, weight="bold"), foreground=ACCENT, spacing1=6, spacing3=8)
            t.tag_configure("h2", font=self._font(size=14, weight="bold"), foreground=INK, spacing1=16, spacing3=5)
            t.tag_configure("h3", font=self._font(size=11, weight="bold"), foreground=ACCENT, spacing1=10, spacing3=3)
            t.tag_configure("bold", font=self._font(size=10, weight="bold"))
            t.tag_configure("italic", font=self._font(size=10, slant="italic"))
            t.tag_configure("bolditalic", font=self._font(size=10, weight="bold", slant="italic"))
            t.tag_configure("code", font=self._font(family=self.MONO, size=9), background="#eef1f5")
            t.tag_configure("codeblock", font=self._font(family=self.MONO, size=9), background="#f1f5f9",
                            lmargin1=16, lmargin2=16, rmargin=16, spacing1=1, spacing3=1, wrap="char")
            t.tag_configure("quote", background="#eaf2ff", lmargin1=18, lmargin2=18, rmargin=18, spacing1=3,
                            spacing3=3, foreground="#1e3a5f")
            t.tag_configure("math", justify="center", spacing1=8, spacing3=10)
            t.tag_configure("muted", foreground=MUTED)
            t.tag_configure("rule", font=self._font(size=2), background="#d1d5db", spacing1=8, spacing3=8)
            for level in (0, 1):
                t.tag_configure(f"list{level}", lmargin1=18 + 22 * level, lmargin2=38 + 22 * level, spacing1=1,
                                spacing3=2)

        # --- navigation ----------------------------------------------------------------------------
        def _on_select(self):
            sel = self.listbox.curselection()
            if sel:
                self.show(PAGES[sel[0]][1])

        def show(self, filename: str, anchor: Optional[str] = None):
            path = os.path.join(self.dir, filename)
            if not os.path.isfile(path):
                self._message(f"Documentation file not found: {path}\n\nThe documentation is on GitHub:\n"
                              f"{GITHUB_DOCS_URL}")
                return
            if filename != self.page:
                with open(path, encoding="utf-8") as fh:
                    self._render(parse_markdown(fh.read(), strict=self.strict))
                self.page = filename
                for idx, (_, f) in enumerate(PAGES):
                    if f == filename:
                        self.listbox.selection_clear(0, "end")
                        self.listbox.selection_set(idx)
            if anchor:
                self.goto(anchor)
            else:
                self.text.yview_moveto(0)

        def goto(self, anchor: str):
            mark = self._anchors.get(anchor)
            if mark:
                self.text.yview(mark)

        def _message(self, msg):
            self.text.configure(state="normal")
            self.text.delete("1.0", "end")
            self.text.insert("1.0", msg)
            self.text.configure(state="disabled")

        def open_link(self, url: str):
            if url.startswith(("http://", "https://")):
                webbrowser.open(url)
            elif url.startswith("#"):
                self.goto(url[1:])
            else:
                name, _, frag = url.partition("#")
                if any(f == name for _, f in PAGES):
                    self.show(name, frag or None)

        # --- rendering ---------------------------------------------------------------------------
        def _render(self, blocks):
            t = self.text
            t.configure(state="normal")
            t.delete("1.0", "end")
            for old in t.mark_names():
                if old.startswith("h-"):
                    t.mark_unset(old)
            self._anchors, self._images = {}, []
            slugs = iter(heading_slugs(blocks))
            for b in blocks:
                if b.kind == "heading":
                    slug = next(slugs)
                    t.mark_set(f"h-{slug}", "end-1c")
                    t.mark_gravity(f"h-{slug}", "left")
                    self._anchors[slug] = f"h-{slug}"
                    self._inline(b.text, [f"h{b.level}"])
                    t.insert("end", "\n")
                elif b.kind == "para":
                    self._inline(b.text, [])
                    t.insert("end", "\n")
                elif b.kind == "list":
                    for level, marker, txt in b.items:
                        t.insert("end", f"{marker}  ", (f"list{level}",))
                        self._inline(txt, [f"list{level}"])
                        t.insert("end", "\n")
                elif b.kind == "quote":
                    t.insert("end", " ", ("quote",))
                    self._inline(b.text, ["quote"])
                    t.insert("end", "\n")
                elif b.kind == "code":
                    t.insert("end", b.text + "\n", ("codeblock",))
                elif b.kind == "math":
                    self._math(b.text, True, ["math"])
                    t.insert("end", "\n", ("math",))
                elif b.kind == "rule":
                    t.insert("end", " \n", ("rule",))
                elif b.kind == "image":
                    self._image(b)
                elif b.kind == "table":
                    self._table(b)
            t.configure(state="disabled")
            t.yview_moveto(0)

        def _math(self, expr, display, tags):
            img = self.math.image(expr, display) if self.math.available else None
            if img is not None:
                self._keep.append(img)
                start = self.text.index("end-1c")
                self.text.image_create("end", image=img, align="center")
                for tag in tags:                  # image_create takes no tags: apply centring / list indent etc.
                    self.text.tag_add(tag, start, "end-1c")
            else:
                self.text.insert("end", latex_to_text(expr), tuple(tags) + ("code",))

        def _inline(self, text, base_tags):
            t = self.text
            for s in parse_inline(text, strict=self.strict):
                style = "bolditalic" if s.bold and s.italic else "bold" if s.bold else "italic" if s.italic else None
                tags = list(base_tags) + ([style] if style else [])
                if s.kind == "math":
                    self._math(s.text, False, tags)
                elif s.kind == "code":
                    t.insert("end", s.text, tuple(tags) + ("code",))
                elif s.kind == "link":
                    name = f"link-{len(self._images)}-{t.index('end')}"
                    self._images.append(name)
                    t.tag_configure(name, foreground="#1f6feb", underline=True)
                    t.tag_bind(name, "<Button-1>", lambda e, u=s.url: self.open_link(u))
                    t.tag_bind(name, "<Enter>", lambda e: t.configure(cursor="hand2"))
                    t.tag_bind(name, "<Leave>", lambda e: t.configure(cursor="arrow"))
                    t.insert("end", s.text, tuple(tags) + (name,))
                else:
                    t.insert("end", s.text, tuple(tags))

        def _image(self, b):
            path = os.path.join(self.dir, b.path)
            try:
                try:
                    from PIL import Image, ImageTk
                    im = Image.open(path)
                    width = max(400, self.text.winfo_width() - 60) if self.text.winfo_width() > 100 else 860
                    if im.width > width:
                        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
                    img = ImageTk.PhotoImage(im)
                except ImportError:
                    img = tk.PhotoImage(file=path)
                    if img.width() > 860:
                        img = img.subsample(-(-img.width() // 860))
                self._keep.append(img)
                self.text.insert("end", "\n", ("math",))
                self.text.image_create("end", image=img)
                self.text.insert("end", "\n", ("math",))
                if b.alt:
                    self.text.insert("end", b.alt + "\n", ("math", "muted"))
            except Exception as e:                                           # missing image: say so, keep going
                self.text.insert("end", f"[image not found: {b.path} ({e})]\n", ("muted",))

        def _table(self, b):
            def clean(s):
                return re.sub(r"[`*]", "", s)
            frame = tk.Frame(self.text, bg="#9aa3af", padx=1, pady=1)
            ncol = len(b.header)
            widths = [max(len(clean(c)) for c in [b.header[k]] + [r[k] for r in b.rows if k < len(r)])
                      for k in range(ncol)]
            total = sum(min(w, 60) for w in widths) or 1
            budget = max(520, self.text.winfo_width() - 80) if self.text.winfo_width() > 100 else 820
            for k in range(ncol):
                px = max(70, int(budget * min(widths[k], 60) / total))
                for r, row in enumerate([b.header] + b.rows):
                    cell = clean(row[k]) if k < len(row) else ""
                    lbl = tk.Label(frame, text=cell, justify="left", anchor="nw", wraplength=px - 12, padx=6, pady=3,
                                   bg="#e5e7eb" if r == 0 else "#ffffff", fg=INK,
                                   font=(self.BODY, 9, "bold" if r == 0 else "normal"))
                    lbl.grid(row=r, column=k, sticky="nsew", padx=(0, 1), pady=(0, 1))
            self._keep.append(frame)
            self.text.insert("end", "\n")
            self.text.window_create("end", window=frame)
            self.text.insert("end", "\n")

    return DocView()
