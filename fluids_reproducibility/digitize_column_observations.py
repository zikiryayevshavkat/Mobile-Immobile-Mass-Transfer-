"""Convert selected Wang et al. (2014) plot markers to numerical observations.

Source DOI: 10.1016/j.crme.2014.07.001. Pixel positions refer to a plot-only
crop from the publisher PDF rendered at 5 pixels per point, resampled to
810 x 360 pixels and placed at (48, 24) on a reading canvas. Crops for
Figures 5a and 6a retain only the lower vertical range to enlarge the markers.
Coordinates are manually selected marker centers, rounded to the nearest pixel.
PDF SHA-256: 480e630798bc24ad9ec3272be6e48ceafd59dd4c66f56ae4b1b603a314913628.
"""

import csv
from pathlib import Path


DESTINATION = Path(__file__).resolve().parent / "results" / "column_observations.csv"
PANELS = (
    # Source, type, feed (g/L), PDF page (zero based), crop, x maximum, full y maximum,
    # displayed y maximum, and marker centers on the reading canvas.
    ("2", "breakthrough", .25, 3, (866, 285, 1997, 897), 3, .20, .20, (
        (48,384),(228,384),(315,314),(360,205),(404,134),(627,60),(715,60),(805,60),
    )),
    ("3a", "breakthrough", .25, 4, (855,286,1985,898), 90, .60, .60, (
        (91,277),(120,277),(149,277),(179,277),(208,273),(223,270),
        (239,252),(253,271),(269,248),(298,243),(327,239),(342,223),(357,216),
        (372,217),(388,198),(418,192),(447,192),(461,187),(476,168),(491,168),
        (506,161),(536,159),(565,154),(580,156),(595,132),(610,131),(625,133),
        (640,132),(655,132),(670,113),(685,113),(700,108),(714,113),(729,109),
        (744,108),(759,108),(773,109),(787,102),
    )),
    ("3b", "deposit", .25, 4, (856,1077,1984,1708), 40, .006, .006, (
        (69,52),(140,154),(241,235),(342,276),(444,304),(545,310),(646,338),(748,338),(828,338),
    )),
    ("5a", "breakthrough", .50, 5, (874,285,1984,912), 90, 1.2, .68, (
        (90,277),(120,265),(150,250),(181,243),(210,236),(225,215),(255,214),
        (285,194),(314,194),(344,182),(374,163),(403,162),(434,141),(463,141),
        (492,135),(522,108),(552,115),(582,108),(611,93),(640,88),(670,83),
        (700,83),(730,78),(759,63),(789,72),
    )),
    ("6a", "deposit", .50, 6, (867,287,1983,926), 40, .030, .010, (
        (67,141),(139,178),(239,216),(341,234),(444,294),(545,305),(645,316),(747,326),(828,341),
    )),
)


def main():
    DESTINATION.parent.mkdir(exist_ok=True)
    with DESTINATION.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("source_figure", "observable", "feed_g_L", "split", "pore_volumes",
                         "distance_m", "observed", "reading_pixel_x", "reading_pixel_y"))
        for panel, kind, feed, _, crop, xmax, ymax, displayed, markers in PANELS:
            height = crop[3]-crop[1]
            displayed = ymax * round(height*displayed/ymax)/height
            for u,v in markers:
                x = (u-48)*xmax/810
                observed = (384-v)*displayed/360
                calibrate = kind=="breakthrough" and feed==.25 and x<=60
                writer.writerow((panel, kind, f"{feed:.2f}", "fit" if calibrate else "holdout",
                                 f"{x:.8f}" if kind=="breakthrough" else "",
                                 f"{x/100:.8f}" if kind=="deposit" else "", f"{observed:.8f}", u,v))


if __name__ == "__main__":
    main()
