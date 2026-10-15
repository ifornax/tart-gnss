"""Campaign table I/O."""

import pandas as pd

from tart_gnss.campaign.output import read_table, write_table


def test_round_trip(tmp_path):
    df = pd.DataFrame({"run_id": ["a"], "stations": [["x", "y"]], "overall_diagnostic_tags": [["NORMAL"]],
                       "evidence_levels": [{"0": 1}], "v": [1.5]})
    p = write_table(df, tmp_path, "t")
    c = read_table(p["csv"])
    assert c["stations"].iloc[0] == ["x", "y"] and c["v"].iloc[0] == 1.5
    if p["parquet"]:
        assert read_table(p["parquet"])["overall_diagnostic_tags"].iloc[0] == ["NORMAL"]
