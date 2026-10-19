"""
Regenerate the screenshots used in the documentation by driving the real window.

    python docs/make_screenshots.py

Windows (it grabs the window rectangle with PIL.ImageGrab) and needs a visible display, numpy and matplotlib. Don't
touch the mouse while it runs: the window is brought to the front and captured. Writes docs/images/calculator.png,
calculator-cu.png, thickogram.png, arxps.png, arxps-fit.png, arxps-fitted.png and depth-trend.png. The inputs are the worked examples quoted in
docs/user-guide.md (the ARXPS ones use the synthetic files in docs/sample/).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from PIL import ImageGrab  # noqa: E402

import app as A  # noqa: E402
import make_sample_data  # noqa: E402
import topofactor as TF  # noqa: E402

IMAGES = os.path.join(HERE, "images")
SIZE = "1500x980+20+20"


def settle(app, n=30):
    for _ in range(n):
        app.update()
        app.after(40)
    app.update_idletasks()
    app.update()


def _window_bounds(app):
    """The window's visible frame (title bar included, no invisible resize border) from the Windows DWM."""
    import ctypes
    from ctypes import wintypes
    hwnd = ctypes.windll.user32.GetParent(app.winfo_id())
    rect = wintypes.RECT()
    DWMWA_EXTENDED_FRAME_BOUNDS = 9
    if ctypes.windll.dwmapi.DwmGetWindowAttribute(hwnd, DWMWA_EXTENDED_FRAME_BOUNDS, ctypes.byref(rect),
                                                 ctypes.sizeof(rect)) != 0:
        raise OSError("DwmGetWindowAttribute failed")
    return rect.left, rect.top, rect.right, rect.bottom


def grab(app, name):
    """Capture exactly the window (so nothing from the desktop behind it ends up in the documentation)."""
    ImageGrab.grab(bbox=_window_bounds(app)).save(os.path.join(IMAGES, name))
    print("wrote", name)


def fill_calculator(app, substrate, overlayers, sub_area, areas):
    s = app.substrate_row
    s.choice_var.set(substrate)
    s._apply_preset(A.SUBSTRATES)
    for i, row in enumerate(app.overlayer_rows):
        key = overlayers[i] if i < len(overlayers) else "None"
        row.choice_var.set(key)
        row._apply_preset(A.OVERLAYERS)
        row.area_var.set(str(areas[i]) if i < len(areas) else "")
    s.area_var.set(str(sub_area))
    app.theta_var.set("90")
    app._on_calculate()


def main():
    os.makedirs(IMAGES, exist_ok=True)
    app = A.App()
    app.geometry(SIZE)
    app.lift()
    app.attributes("-topmost", True)
    settle(app, 20)
    app._place_sash()
    settle(app, 10)
    nb = app.nametowidget(app.calc_tab.winfo_parent())

    nb.select(app.calc_tab)
    fill_calculator(app, "Al (metal)", ["Al2O3"], 35, [65])
    settle(app)
    grab(app, "calculator.png")

    fill_calculator(app, "Cu (metal, 2p3/2)", ["CuO (2p3/2)", "Cu2O (2p3/2)", "Cu(OH)2 (2p3/2)"], 25, [30, 40, 13])
    settle(app)
    grab(app, "calculator-cu.png")

    nb.select(app.thick_tab)
    tab = app.thickogram
    for key, value in {"Eo": "285", "Es": "99.5", "Io": "800", "Is": "3000", "So": "1.0", "Ss": "3.0",
                       "Lambda": "3.0", "Theta": "0", "Esource": "1486.6"}.items():
        tab.vars[key].set(value)
    tab.shape_var.set(TF.SHAPES["cylinder"].label)
    tab.calculate()
    settle(app)
    grab(app, "thickogram.png")

    nb.select(app.arxps_tab)                       # the ARXPS worked example: the synthetic sample file
    tab = app.arxps
    tab.load_file(make_sample_data.SAMPLE)
    tab.v["region"].set("Si2p Scan")
    tab._region_changed()
    tab.sub_preset.set("Si (metal/elemental)")
    tab._apply_preset("sub")
    tab.ov_preset.set("SiO2")
    tab._apply_preset("ov")
    tab.auto_windows()
    tab.calculate()
    settle(app)
    tab.plot_nb.select(0)
    settle(app)
    grab(app, "arxps.png")
    tab.plot_nb.select(1)
    settle(app)
    grab(app, "arxps-fit.png")

    tab.load_file(make_sample_data.SAMPLE_FITTED)     # the same series with a CasaXPS-style fit stored in it
    tab.v["region"].set("Si2p Scan")
    tab._region_changed()
    tab.sub_preset.set("Si (metal/elemental)")
    tab._apply_preset("sub")
    tab.ov_preset.set("SiO2")
    tab._apply_preset("ov")
    tab.calculate()
    tab.setup_nb.select(2)                            # the Components page
    tab.plot_nb.select(0)
    settle(app)
    grab(app, "arxps-fitted.png")

    nb.select(app.depth_tab)                          # the Depth profile tab on the fitted synthetic series
    tab = app.depth
    tab.load_file(make_sample_data.SAMPLE_FITTED)
    settle(app)
    grab(app, "depth-trend.png")
    app.destroy()


if __name__ == "__main__":
    main()
