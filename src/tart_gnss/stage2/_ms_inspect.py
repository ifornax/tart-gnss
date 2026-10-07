"""Summarise a CASA Measurement Set as JSON (standalone; needs python-casacore).

Kept free of tart_gnss imports so it can run under another interpreter
(e.g. the one that provides tart2ms):  ``python _ms_inspect.py <ms>``.
"""

import json
import sys

MJD_UNIX_OFFSET_S = 3506716800.0  # seconds from MJD 0 (1858-11-17) to 1970-01-01


def summarise(path):
    """Return a dict describing the MS at ``path`` (raises on unreadable tables)."""
    from casacore.tables import table

    out = {"path": path}
    with table(path, ack=False) as main:
        cols = main.colnames()
        out["columns"] = cols
        out["number_of_rows"] = int(main.nrows())
        for name in ("TIME", "UVW", "DATA", "FLAG"):
            out["has_" + name.lower()] = name in cols
        if out["number_of_rows"] and "TIME" in cols:
            t = main.getcol("TIME")
            out["first_time_unix"] = float(t.min()) - MJD_UNIX_OFFSET_S
            out["last_time_unix"] = float(t.max()) - MJD_UNIX_OFFSET_S
            out["number_of_times"] = int(len(set(t.tolist())))
        if out["number_of_rows"] and "FLAG" in cols:
            out["flagged_fraction"] = float(main.getcol("FLAG").mean())
        if out["number_of_rows"] and "UVW" in cols:
            uvw = main.getcol("UVW")
            out["max_baseline_m"] = float((uvw ** 2).sum(axis=1).max() ** 0.5)
    keywords = set(table(path, ack=False).getkeywords())
    for sub in ("ANTENNA", "SPECTRAL_WINDOW", "FIELD", "OBSERVATION"):
        out["has_" + sub.lower() + "_table"] = sub in keywords
    if out["has_antenna_table"]:
        with table(path + "/ANTENNA", ack=False) as ant:
            out["number_of_antennas"] = int(ant.nrows())
    if out["has_spectral_window_table"]:
        with table(path + "/SPECTRAL_WINDOW", ack=False) as spw:
            out["number_of_spectral_windows"] = int(spw.nrows())
            out["number_of_channels"] = [int(x) for x in spw.getcol("NUM_CHAN")]
            out["channel_frequencies_hz"] = [float(x) for x in spw.getcol("CHAN_FREQ").ravel()]
            out["total_bandwidth_hz"] = [float(x) for x in spw.getcol("TOTAL_BANDWIDTH")]
    if out["has_observation_table"]:
        with table(path + "/OBSERVATION", ack=False) as obs:
            out["telescope_name"] = [str(x) for x in obs.getcol("TELESCOPE_NAME")]
    if out["has_field_table"]:
        with table(path + "/FIELD", ack=False) as field:
            out["number_of_fields"] = int(field.nrows())
    return out


if __name__ == "__main__":
    print(json.dumps(summarise(sys.argv[1])))
