"""
Documentation tests: docs/*.md must stay within the Markdown subset the in-app help can show, every equation must be
typesettable, links/images must resolve, citations must exist, and the worked examples must match the code.

    python -m unittest discover -s tests -v
"""
import glob
import importlib.util
import os
import re
import sys
import unittest

ROOT = os.path.join(os.path.dirname(__file__), os.pardir)
sys.path.insert(0, ROOT)

import app as A  # noqa: E402
import arxps as AR  # noqa: E402
import help_viewer as H  # noqa: E402
import thickogram as T  # noqa: E402
import topofactor as TFAC  # noqa: E402
import vamas as VM  # noqa: E402

DOCS = os.path.join(ROOT, "docs")

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.mathtext import MathTextParser
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


def _read(name):
    with open(os.path.join(DOCS, name), encoding="utf-8") as fh:
        return fh.read()


def _all_pages():
    return [os.path.basename(p) for p in sorted(glob.glob(os.path.join(DOCS, "*.md")))]


class Structure(unittest.TestCase):
    def test_the_app_pages_exist_and_are_the_docs_folder_pages(self):
        names = [f for _, f in H.PAGES]
        for f in names:
            self.assertTrue(os.path.isfile(os.path.join(DOCS, f)), f)
        self.assertEqual(sorted(names), _all_pages())

    def test_every_page_parses_in_strict_mode(self):
        for f in _all_pages():
            blocks = H.parse_markdown(_read(f), strict=True)
            self.assertTrue(blocks and blocks[0].kind == "heading" and blocks[0].level == 1, f)
            for text in H.iter_inline_texts(blocks):
                H.parse_inline(text, strict=True)             # raises on anything unsupported

    def test_images_exist_and_are_big_enough(self):
        for f in _all_pages():
            for b in H.parse_markdown(_read(f)):
                if b.kind == "image":
                    path = os.path.join(DOCS, b.path)
                    self.assertTrue(os.path.isfile(path), f"{f}: missing image {b.path}")
                    with open(path, "rb") as fh:
                        head = fh.read(24)
                    self.assertEqual(head[:8], b"\x89PNG\r\n\x1a\n", f"{b.path} is not a PNG")
                    width = int.from_bytes(head[16:20], "big")
                    self.assertGreaterEqual(width, 700, f"{b.path} is only {width}px wide")
                    self.assertTrue(b.alt, f"{b.path} needs alt text")

    def test_links_resolve_including_anchors(self):
        slugs = {f: H.heading_slugs(H.parse_markdown(_read(f))) for f in _all_pages()}
        checked = 0
        for f in _all_pages():
            blocks = H.parse_markdown(_read(f))
            for text in H.iter_inline_texts(blocks):
                for span in H.parse_inline(text):
                    if span.kind != "link" or span.url.startswith(("http://", "https://")):
                        continue
                    name, _, frag = span.url.partition("#")
                    target = name or f
                    self.assertIn(target, slugs, f"{f}: link to unknown page {span.url}")
                    if frag:
                        self.assertIn(frag, slugs[target], f"{f}: no heading #{frag} in {target}")
                    checked += 1
        self.assertGreater(checked, 10)


@unittest.skipUnless(HAVE_MPL, "matplotlib not installed")
class Equations(unittest.TestCase):
    def test_every_equation_is_typesettable_and_github_safe(self):
        parser = MathTextParser("agg")
        forbidden = (r"\begin", r"\tfrac", r"\newcommand", r"\\\\")
        count = 0
        for f in _all_pages():
            for expr in H.iter_math(H.parse_markdown(_read(f))):
                self.assertTrue(expr.strip(), f"{f}: empty equation")
                for bad in forbidden:
                    self.assertNotIn(bad, expr, f"{f}: {bad} is not supported by both GitHub and the app: {expr}")
                try:
                    parser.parse("$" + expr + "$", dpi=100)
                except Exception as e:                                              # noqa: BLE001
                    self.fail(f"{f}: mathtext cannot typeset {expr!r}: {e}")
                count += 1
        self.assertGreater(count, 30)


class Citations(unittest.TestCase):
    def test_reference_list_has_the_papers_with_verified_dois(self):
        refs = _read("references.md")
        for doi in ("10.1002/sia.3044",                    # Shard, Wang, Spencer 2009 (Topofactors)
                    "10.1002/sia.740150109",               # Strohmeier 1990
                    "10.1002/sia.6239",                    # Biesinger 2017
                    "10.1002/sia.740210302",               # Tanuma, Powell, Penn 1994 (TPP-2M)
                    "10.1016/0368-2048(94)02270-4",        # Cumpson 1995 (angle-resolved XPS)
                    "10.1002/(SICI)1096-9918(199706)25:6<430::AID-SIA254>3.0.CO;2-7",   # Cumpson, Seah 1997
                    "10.1103/PhysRevB.5.4709",             # Shirley 1972
                    "10.1002/sia.740130202",               # Dench, Hazell, Seah 1988 (VAMAS)
                    "10.1002/1096-9918(200006)29:6<403::AID-SIA884>3.0.CO;2-8"):   # Cumpson 2000
            self.assertTrue(doi in refs, f"references.md lacks the DOI {doi}")
        self.assertTrue("XPS topofactors: determining overlayer thickness on particles and" in refs)
        self.assertTrue("Shard, J. Wang and S. J. Spencer" in refs)

    def test_numbered_citations_point_at_real_references(self):
        refs = _read("references.md")
        numbers = {int(m) for m in re.findall(r"(?m)^(\d+)\. ", refs)}
        self.assertEqual(numbers, set(range(1, max(numbers) + 1)))
        for f in _all_pages():
            blocks = H.parse_markdown(_read(f))
            for text in H.iter_inline_texts(blocks):
                plain = "".join(s.text for s in H.parse_inline(text) if s.kind == "text")
                for n in re.findall(r"\[(\d+)\]", plain):
                    self.assertIn(int(n), numbers, f"{f}: citation [{n}] has no entry in references.md")

    def test_the_topofactor_paper_is_cited_where_the_feature_is_documented(self):
        for f in ("methods.md", "user-guide.md", "references.md"):
            self.assertTrue(re.search(r"\[4\]|Shard", _read(f)), f"{f} does not cite the Topofactor paper")
        for name in ("NOTICE", "README.md"):
            with open(os.path.join(ROOT, name), encoding="utf-8") as fh:
                text = fh.read()
            self.assertTrue("Shard" in text and "10.1002/sia.3044" in text, f"{name} must cite Shard et al. 2009")


class ClaimsMatchTheCode(unittest.TestCase):
    """Numbers quoted in the docs are recomputed here, so an edit to the maths can't leave the docs stale."""

    def test_al_example(self):
        sub, ov = A.SUBSTRATES["Al (metal)"], A.OVERLAYERS["Al2O3"]
        n = lambda c: A.volume_density(c.density_g_cm3, c.molar_mass_g_mol, c.atoms_per_formula)   # noqa: E731
        res = A.strohmeier_thickness(90, n(sub), sub.imfp_nm, 1.0, 35, [A.OverlayerInput("x", n(ov), ov.imfp_nm, 1.0, 65)])
        self.assertIn(f"{res['d_nm']:.2f} nm ({res['d_angstrom']:.1f} angstrom)", _read("user-guide.md"))

    def test_cu_lmm_example(self):
        sub = A.SUBSTRATES["Cu (metal, LMM)"]
        comps = []
        for key, area in (("Cu2O (LMM)", 63), ("Cu(OH)2 (LMM)", 13)):
            o = A.OVERLAYERS[key]
            comps.append(A.OverlayerInput(key, A.volume_density(o.density_g_cm3, o.molar_mass_g_mol, o.atoms_per_formula),
                                          o.imfp_nm, 1.0, area))
        res = A.strohmeier_thickness(90, A.volume_density(sub.density_g_cm3, sub.molar_mass_g_mol, 1), sub.imfp_nm, 1.0,
                                     25, comps)
        self.assertIn(f"**{res['d_nm']:.2f} nm**", _read("user-guide.md"))

    def test_topofactor_example(self):
        r = T.thickness(io=800.0, is_=3000.0, so=1.0, ss=3.0, be_o=285.0, be_s=99.5, source_ev=1486.6, lam_nm=3.0,
                        theta_normal_deg=0.0)
        c = TFAC.apply_topofactor(r, "cylinder", "eq6")
        guide = _read("user-guide.md")
        self.assertIn(f"**{r.t_nm:.2f} nm**", guide)
        self.assertIn(f"**{c.tf:.3f}**", guide)
        self.assertIn(f"**{c.t_nm:.2f} nm**", guide)

    def test_arxps_worked_example_and_sample_file(self):
        sample = os.path.join(DOCS, "sample", "si-native-oxide-arxps.vms")
        sub, ov = A.SUBSTRATES["Si (metal/elemental)"], A.OVERLAYERS["SiO2"]
        n = lambda c: A.volume_density(c.density_g_cm3, c.molar_mass_g_mol, c.atoms_per_formula)   # noqa: E731
        rinf = ov.rsf * n(ov) * ov.imfp_nm / (sub.rsf * n(sub) * sub.imfp_nm)
        series = VM.read(sample).series("Si2p Scan", "exp1")[None]
        w = AR.auto_two_state_windows([(b.be_ev, AR.intensity(b)) for _, b in series], sub.be_ev, ov.be_ev)
        res = AR.analyse_two_states(series, w, ov.imfp_nm, rinf)
        guide = _read("user-guide.md")
        self.assertIn(f"**{res.d_nm:.2f} nm**", guide)
        self.assertIn(f"the envelope is {w.env_lo:.1f} to {w.env_hi:.1f} eV, split at {w.split:.1f} eV", guide)
        self.assertIn(f"lambda {ov.imfp_nm:g} nm, $R_\\infty$ = {rinf:.3f}", guide)
        self.assertAlmostEqual(res.d_nm, 1.5, delta=0.02)               # "within 0.02 nm of the value put in"
        self.assertEqual([r.tilt for r in res.rows if not r.used], [75.0])
        self.assertAlmostEqual(res.line[0], 0.0, delta=2 * res.line[2])

    def test_the_sample_file_is_what_the_generator_writes(self):
        import tempfile
        spec = importlib.util.spec_from_file_location("make_sample_data", os.path.join(DOCS, "make_sample_data.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as d:
            fresh = os.path.join(d, "s.vms")
            import synthetic
            synthetic.write_sample(fresh, d_nm=mod.TRUE_THICKNESS_NM)
            with open(fresh, "rb") as a, open(mod.SAMPLE, "rb") as b:
                self.assertEqual(a.read(), b.read(), "docs/sample is stale: run  python docs/make_sample_data.py")

    def test_table_1_and_the_constants_in_methods(self):
        methods = _read("methods.md")
        labels = {"cylinder": "Cylinder (fibre)", "sphere": "Sphere, X-rays reach all of it",
                  "sphere_shadowed": "Sphere, X-ray shadowed"}
        rows = {}
        for line in methods.splitlines():
            if line.startswith("| ") and line.count("|") >= 10:
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                rows[cells[0]] = cells[1:]
        for key, label in labels.items():
            self.assertIn(label, rows, f"methods.md has no Table 1 row for {label}")
            want = list(TFAC.TABLE_1[key]) + [round(TFAC.SHAPES[key].thin_limit, 2), TFAC.SHAPES[key].constant]
            got = [float(c) for c in rows[label]]
            self.assertEqual(got, want, f"methods.md Table 1 row for {label} disagrees with topofactor.py")

    def test_the_quoted_validity_numbers(self):
        methods = _read("methods.md")
        self.assertIn(f"$0.1 < t < {TFAC.T_MAX:g}$", methods)
        self.assertIn(f"between {TFAC.B_USUAL[0]} and {TFAC.B_USUAL[1]}", methods)

    def test_presets_page_is_up_to_date(self):
        spec = importlib.util.spec_from_file_location("make_presets", os.path.join(DOCS, "make_presets.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        self.assertEqual(_read("presets.md"), mod.build(),
                         "docs/presets.md is stale: run  python docs/make_presets.py")


class MarkdownParser(unittest.TestCase):
    def test_block_types(self):
        md = ("# T\n\npara one\ncontinues\n\n- a\n- b\n  more\n  - nested\n\n1. x\n2. y\n\n> quote\n\n```\ncode\n```\n\n"
              "```math\nx^2\n```\n\n| a | b |\n|---|---|\n| 1 | 2 |\n\n![alt](images/x.png)\n\n---\n")
        kinds = [b.kind for b in H.parse_markdown(md, strict=True)]
        self.assertEqual(kinds, ["heading", "para", "list", "list", "quote", "code", "math", "table", "image", "rule"])
        blocks = H.parse_markdown(md)
        self.assertEqual(blocks[1].text, "para one continues")
        self.assertEqual([(l, m) for l, m, _ in blocks[2].items], [(0, "•"), (0, "•"), (1, "•")])
        self.assertIn("b more", blocks[2].items[1][2])
        self.assertEqual(blocks[7].header, ["a", "b"])
        self.assertEqual(blocks[7].rows, [["1", "2"]])

    def test_inline_spans(self):
        spans = H.parse_inline("a **bold *both* x** `c` $x_1$ [l](u.md#h) https://e.org/z. end", strict=True)
        kinds = [(s.kind, s.text, s.bold, s.italic) for s in spans]
        self.assertIn(("math", "x_1", False, False), kinds)
        self.assertIn(("code", "c", False, False), kinds)
        self.assertIn(("text", "both", True, True), kinds)
        self.assertIn(("link", "l", False, False), kinds)
        link = [s for s in spans if s.kind == "link" and s.text.startswith("https")][0]
        self.assertEqual(link.url, "https://e.org/z")                       # trailing full stop is not part of the URL

    def test_escapes_and_literal_asterisks(self):
        spans = H.parse_inline(r"2 \* 3 and \$5", strict=True)
        self.assertEqual("".join(s.text for s in spans), "2 * 3 and $5")

    def test_strict_mode_rejects_unsupported_markdown(self):
        for bad in ("<b>x</b>", "~~x~~", "an _emphasis_ here", "unclosed **bold", "price $5 only", "`unclosed"):
            with self.assertRaises(ValueError, msg=bad):
                H.parse_inline(bad, strict=True)
        for bad in ("#### too deep", "```\nunclosed", "text\n---\n", "| a | b |\n|---|---|\n| 1 |\n",
                    "| a |\n|---|\n| $x$ |\n", "- a\n    - b\n        - c\n"):
            with self.assertRaises(ValueError, msg=bad):
                H.parse_markdown(bad, strict=True)

    def test_anchors_follow_github(self):
        self.assertEqual(H.slugify("Strohmeier equation (Calculator tab)"), "strohmeier-equation-calculator-tab")
        self.assertEqual(H.slugify("Which tab should I use?"), "which-tab-should-i-use")
        self.assertEqual(H.slugify("IMFPs from TPP-2M"), "imfps-from-tpp-2m")
        blocks = H.parse_markdown("# A\n\n## B\n\n## B\n")
        self.assertEqual(H.heading_slugs(blocks), ["a", "b", "b-1"])


class PlainTextEquations(unittest.TestCase):
    """The no-matplotlib fallback must at least be readable."""

    def test_samples(self):
        self.assertEqual(H.latex_to_text(r"\frac{a}{b}"), "(a)/(b)")
        self.assertEqual(H.latex_to_text(r"\lambda_{ox}\,\sin\theta"), "λ_ox sinθ")
        self.assertEqual(H.latex_to_text(r"x^{2}"), "x²")
        self.assertEqual(H.latex_to_text(r"\frac{2}{\pi}\approx 0.64"), "(2)/(π)≈ 0.64".replace("≈ ", "≈ "))
        self.assertEqual(H.latex_to_text(r"\left(N_m/N_{ox}\right)"), "(N_m/N_ox)")
        self.assertNotIn("\\", H.latex_to_text(r"\mathrm{TF}\approx c_0+\left(c_1+c_2\ln B\right)t^{k}"))

    def test_every_equation_in_the_docs_converts_without_backslashes(self):
        for f in _all_pages():
            for expr in H.iter_math(H.parse_markdown(_read(f))):
                text = H.latex_to_text(expr)
                self.assertTrue(text, expr)
                self.assertNotIn("\\", text, expr)


if __name__ == "__main__":
    unittest.main()
