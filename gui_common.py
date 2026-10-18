"""
Shared GUI helpers for the XPS Thickness Suite: the palette, the ttk theme/styles (`setup_style`), the hover
`Tooltip` and `_to_float`. Kept in its own module so that app.py and the tab modules (e.g. thickogram_tab.py)
can both use them without importing each other.
"""
import tkinter as tk
from tkinter import ttk

# Palette (one place, so the look is easy to tweak)
BG = "#eef1f5"
CARD = "#ffffff"
INK = "#1f2937"
MUTED = "#6b7280"
ACCENT = "#1f6feb"
ACCENT_HOVER = "#1a5fd0"
ACCENT_PRESSED = "#174fb0"
BORDER = "#8d96a3"
BTN = "#e5e7eb"
BTN_HOVER = "#d6dbe2"
BTN_PRESSED = "#c3cad3"
BADGE = {   # kind -> (background, foreground)
    "idle": ("#e5e7eb", "#4b5563"),
    "ok": ("#d1fae5", "#065f46"),
    "warn": ("#fef3c7", "#92400e"),
    "bad": ("#fee2e2", "#991b1b"),
}


def setup_style(root, compact=False):
    """Switch to the 'clam' theme and define the button/card/entry styles used by the app."""
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    base = 8 if compact else 9
    font = ("Segoe UI", base)
    root.configure(bg=BG)

    style.configure(".", background=BG, foreground=INK, font=font, bordercolor=BORDER, focuscolor=ACCENT)
    style.configure("TFrame", background=BG)
    style.configure("Card.TFrame", background=CARD)
    style.configure("TLabel", background=BG, foreground=INK)
    style.configure("Card.TLabel", background=CARD, foreground=INK)
    style.configure("Muted.TLabel", background=BG, foreground=MUTED)
    style.configure("CardMuted.TLabel", background=CARD, foreground=MUTED)
    style.configure("Header.TLabel", background=BG, foreground=INK, font=("Segoe UI", base + 4, "bold"))
    style.configure("ColHead.TLabel", background=CARD, foreground=MUTED, font=("Segoe UI", base - 1, "bold"))
    style.configure("RowHead.TLabel", background=CARD, foreground=INK, font=("Segoe UI", base, "bold"))
    style.configure("Card.TCheckbutton", background=CARD, foreground=INK)
    style.map("Card.TCheckbutton", background=[("active", CARD)])

    style.configure("Card.TLabelframe", background=CARD, bordercolor=BORDER, relief="solid", borderwidth=1)
    style.configure("Card.TLabelframe.Label", background=CARD, foreground=ACCENT, font=("Segoe UI", base + 1, "bold"))
    style.configure("TSeparator", background="#d1d5db")

    # Buttons: a clearly bordered, raised secondary style and a bold blue primary style, with hover/pressed states.
    style.configure("TButton", background=BTN, foreground=INK, bordercolor=BORDER, lightcolor="#ffffff",
                    darkcolor="#c1c7cf", relief="raised", borderwidth=1, padding=(10, 5), focuscolor=BTN, font=font)
    style.map("TButton",
              background=[("disabled", "#f3f4f6"), ("pressed", BTN_PRESSED), ("active", BTN_HOVER)],
              foreground=[("disabled", "#9ca3af")],
              relief=[("pressed", "sunken"), ("!pressed", "raised")])
    style.configure("Primary.TButton", background=ACCENT, foreground="#ffffff", bordercolor=ACCENT_PRESSED,
                    lightcolor="#4a8df0", darkcolor=ACCENT_PRESSED, relief="raised", borderwidth=1,
                    padding=(18, 9), focuscolor=ACCENT, font=("Segoe UI", base + 2, "bold"))
    style.map("Primary.TButton",
              background=[("disabled", "#9db9e8"), ("pressed", ACCENT_PRESSED), ("active", ACCENT_HOVER)],
              foreground=[("disabled", "#ffffff")],
              relief=[("pressed", "sunken"), ("!pressed", "raised")])

    # Input fields
    style.configure("TEntry", fieldbackground="#ffffff", bordercolor=BORDER, lightcolor=BORDER, darkcolor=BORDER,
                    padding=3)
    style.map("TEntry", bordercolor=[("focus", ACCENT)], lightcolor=[("focus", ACCENT)],
              darkcolor=[("focus", ACCENT)])
    style.configure("TCombobox", fieldbackground="#ffffff", bordercolor=BORDER, lightcolor=BORDER, darkcolor=BORDER,
                    arrowcolor=INK, padding=3)
    style.map("TCombobox", fieldbackground=[("readonly", "#ffffff")], selectbackground=[("readonly", "#ffffff")],
              selectforeground=[("readonly", INK)], bordercolor=[("focus", ACCENT)])

    style.configure("TNotebook", background=BG, borderwidth=0, tabmargins=(2, 4, 2, 0))
    style.configure("TNotebook.Tab", background="#dfe4ea", foreground=MUTED, padding=(18, 7),
                    font=("Segoe UI", base + 1, "bold"))
    style.map("TNotebook.Tab", background=[("selected", CARD), ("active", "#e9edf2")],
              foreground=[("selected", ACCENT)], padding=[("selected", (18, 9))])

    style.configure("Horizontal.TScale", background=CARD, troughcolor="#d1d5db")
    style.configure("TPanedwindow", background=BG)
    style.configure("Sash", sashthickness=10, gripcount=16, handlesize=14)


class Tooltip:
    """Small hover tooltip: appears after a short delay, disappears on leave/click."""

    def __init__(self, widget, text, delay=450):
        self.widget, self.text, self.delay = widget, text, delay
        self._tip = None
        self._after = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._after = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after is not None:
            self.widget.after_cancel(self._after)
            self._after = None

    def _show(self):
        if self._tip is not None or not self.text:
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self._tip, text=self.text, justify="left", background="#111827", foreground="#ffffff",
                 padx=8, pady=4, font=("Segoe UI", 8), wraplength=320).pack()

    def _hide(self, _event=None):
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None


def _to_float(text):
    """Float from an entry's text, or None if blank/invalid."""
    try:
        return float(str(text).strip())
    except ValueError:
        return None
