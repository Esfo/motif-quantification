#!/usr/bin/env python3

import argparse
import json
import math
import sqlite3
from pathlib import Path


DEFAULT_LIMIT = 20


def connect(path):
    path = Path(path)

    if not path.exists():
        raise SystemExit(f"missing sqlite file: {path}")

    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    return conn


def table_exists(conn, table):
    row = conn.execute(
        """
        SELECT 1
        FROM sqlite_master
        WHERE type = 'table'
          AND name = ?
        """,
        (table,),
    ).fetchone()

    return row is not None


def column_exists(conn, table, column):
    if not table_exists(conn, table):
        return False

    rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
    return column in {row["name"] for row in rows}


def count_rows(conn, table):
    if not table_exists(conn, table):
        return None

    return conn.execute(f"SELECT count(*) AS n FROM {table}").fetchone()["n"]


def fetch_parameters(conn):
    if not table_exists(conn, "parameters"):
        return {}

    params = {}

    for row in conn.execute("SELECT key, value FROM parameters ORDER BY key"):
        key = row["key"]
        value = row["value"]

        try:
            params[key] = json.loads(value)
        except Exception:
            params[key] = value

    return params


def quantiles(values, probs=(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0)):
    values = sorted(v for v in values if v is not None and math.isfinite(v))

    if not values:
        return {}

    n = len(values)
    out = {}

    for p in probs:
        if n == 1:
            value = values[0]
        else:
            pos = p * (n - 1)
            left = int(math.floor(pos))
            right = int(math.ceil(pos))

            if left == right:
                value = values[left]
            else:
                frac = pos - left
                value = values[left] * (1.0 - frac) + values[right] * frac

        label = f"p{int(p * 100):02d}" if p not in (0, 1.0) else ("min" if p == 0 else "max")
        out[label] = value

    return out


def fetch_values(conn, sql):
    return [row[0] for row in conn.execute(sql)]


def print_header(title):
    print()
    print(title)
    print("=" * len(title))


def format_value(value):
    if value is None:
        return ""

    if isinstance(value, int):
        return str(value)

    if isinstance(value, float):
        if abs(value) >= 1000:
            return f"{value:.3f}"
        if abs(value) >= 1:
            return f"{value:.6g}"
        return f"{value:.6g}"

    return str(value)


def print_table(headers, rows):
    rows = list(rows)

    if not rows:
        print("(none)")
        return

    str_rows = [[format_value(value) for value in row] for row in rows]
    str_headers = [str(header) for header in headers]

    widths = [
        max(len(str_headers[i]), *(len(row[i]) for row in str_rows))
        for i in range(len(str_headers))
    ]

    print("  ".join(str_headers[i].ljust(widths[i]) for i in range(len(str_headers))))
    print("  ".join("-" * widths[i] for i in range(len(str_headers))))

    for row in str_rows:
        print("  ".join(row[i].rjust(widths[i]) for i in range(len(str_headers))))


def print_key_values(rows):
    rows = [(str(k), format_value(v)) for k, v in rows]
    width = max((len(k) for k, _ in rows), default=0)

    for key, value in rows:
        print(f"{key.ljust(width)} : {value}")


def print_quantile_block(title, values, unit=""):
    q = quantiles(values)

    print_header(title)

    if not q:
        print("(none)")
        return

    rows = []

    for key in ("min", "p01", "p05", "p25", "p50", "p75", "p95", "p99", "max"):
        value = q.get(key)

        if value is not None and unit:
            rows.append((key, f"{format_value(value)} {unit}"))
        else:
            rows.append((key, value))

    print_key_values(rows)


def table_counts(conn):
    tables = [
        "parameters",
        "scans",
        "lines",
        "features",
        "feature_traces",
        "isotope_edges",
        "distributions",
        "distribution_members",
        "analytes",
        "analyte_members",
    ]

    rows = []

    for table in tables:
        n = count_rows(conn, table)

        if n is not None:
            rows.append((table, n))

    return rows


def metadata_rows(db_path, params):
    path = Path(db_path)

    rows = [
        ("sqlite", str(path)),
        ("size_mb", path.stat().st_size / 1_000_000),
    ]

    for key in ("script", "script_version", "mzml", "store_edges", "workers"):
        if key in params:
            rows.append((key, params[key]))

    if "counts" in params and isinstance(params["counts"], dict):
        counts = params["counts"]

        for key in sorted(counts):
            rows.append((f"stored_count.{key}", counts[key]))

    return rows


def print_distribution_by_members(conn):
    print_header("Distributions by isotope member count")

    rows = conn.execute(
        """
        SELECT n_members, count(*) AS n
        FROM distributions
        GROUP BY n_members
        ORDER BY n_members
        """
    ).fetchall()

    print_table(
        ["n_members", "count", "fraction"],
        [
            (
                row["n_members"],
                row["n"],
                row["n"] / max(1, count_rows(conn, "distributions")),
            )
            for row in rows
        ],
    )


def print_distribution_by_charge(conn):
    print_header("Distributions by charge")

    rows = conn.execute(
        """
        SELECT charge, count(*) AS n
        FROM distributions
        GROUP BY charge
        ORDER BY charge
        """
    ).fetchall()

    total = max(1, count_rows(conn, "distributions"))

    print_table(
        ["charge", "count", "fraction"],
        [(row["charge"], row["n"], row["n"] / total) for row in rows],
    )


def print_distribution_charge_members(conn):
    print_header("Distributions by charge and isotope member count")

    rows = conn.execute(
        """
        SELECT charge, n_members, count(*) AS n
        FROM distributions
        GROUP BY charge, n_members
        ORDER BY charge, n_members
        """
    ).fetchall()

    print_table(
        ["charge", "n_members", "count"],
        [(row["charge"], row["n_members"], row["n"]) for row in rows],
    )


def print_distribution_score_summaries(conn):
    print_quantile_block(
        "Distribution score quantiles",
        fetch_values(conn, "SELECT score FROM distributions"),
    )

    print_quantile_block(
        "Distribution quality quantiles",
        fetch_values(conn, "SELECT quality FROM distributions"),
    )

    print_quantile_block(
        "Distribution RT width quantiles",
        fetch_values(conn, "SELECT rt_end - rt_start FROM distributions"),
        unit="min",
    )

    print_quantile_block(
        "Distribution neutral mass quantiles",
        fetch_values(conn, "SELECT neutral_mass FROM distributions"),
        unit="Da",
    )


def print_feature_summaries(conn):
    if not table_exists(conn, "features"):
        return

    print_quantile_block(
        "Feature point-count quantiles",
        fetch_values(conn, "SELECT n_points FROM features"),
    )

    print_quantile_block(
        "Feature RT width quantiles",
        fetch_values(conn, "SELECT rt_end - rt_start FROM features"),
        unit="min",
    )

    print_quantile_block(
        "Feature height quantiles",
        fetch_values(conn, "SELECT height FROM features"),
    )

    print_quantile_block(
        "Feature area quantiles",
        fetch_values(conn, "SELECT area FROM features"),
    )


def print_analyte_summaries(conn):
    if not table_exists(conn, "analytes"):
        return

    print_quantile_block(
        "Analyte distribution-count quantiles",
        fetch_values(conn, "SELECT n_distributions FROM analytes"),
    )

    print_quantile_block(
        "Analyte RT width quantiles",
        fetch_values(conn, "SELECT rt_end - rt_start FROM analytes"),
        unit="min",
    )

    print_quantile_block(
        "Analyte neutral mass quantiles",
        fetch_values(conn, "SELECT neutral_mass FROM analytes"),
        unit="Da",
    )

    print_header("Analytes by charge span")

    rows = conn.execute(
        """
        SELECT charge_min, charge_max, count(*) AS n
        FROM analytes
        GROUP BY charge_min, charge_max
        ORDER BY charge_min, charge_max
        """
    ).fetchall()

    print_table(
        ["charge_min", "charge_max", "count"],
        [(row["charge_min"], row["charge_max"], row["n"]) for row in rows],
    )


def suspicious_tiny_rt(conn, limit, max_width):
    print_header(f"Suspicious distributions: tiny RT width <= {max_width:g} min")

    rows = conn.execute(
        """
        SELECT
            distribution_id,
            charge,
            n_members,
            neutral_mass,
            mono_mz,
            rt_start,
            rt_apex,
            rt_end,
            rt_end - rt_start AS rt_width,
            score,
            quality,
            ms1_start,
            ms1_apex,
            ms1_end
        FROM distributions
        WHERE rt_end - rt_start <= ?
        ORDER BY quality ASC, score ASC
        LIMIT ?
        """,
        (max_width, limit),
    ).fetchall()

    print_table(
        [
            "distribution_id",
            "z",
            "members",
            "mass",
            "mono_mz",
            "rt_width",
            "score",
            "quality",
            "ms1",
        ],
        [
            (
                row["distribution_id"],
                row["charge"],
                row["n_members"],
                row["neutral_mass"],
                row["mono_mz"],
                row["rt_width"],
                row["score"],
                row["quality"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def suspicious_huge_rt(conn, limit, min_width):
    print_header(f"Suspicious distributions: huge RT width >= {min_width:g} min")

    rows = conn.execute(
        """
        SELECT
            distribution_id,
            charge,
            n_members,
            neutral_mass,
            mono_mz,
            rt_start,
            rt_apex,
            rt_end,
            rt_end - rt_start AS rt_width,
            score,
            quality,
            ms1_start,
            ms1_apex,
            ms1_end
        FROM distributions
        WHERE rt_end - rt_start >= ?
        ORDER BY rt_width DESC, quality ASC
        LIMIT ?
        """,
        (min_width, limit),
    ).fetchall()

    print_table(
        [
            "distribution_id",
            "z",
            "members",
            "mass",
            "mono_mz",
            "rt_width",
            "score",
            "quality",
            "ms1",
        ],
        [
            (
                row["distribution_id"],
                row["charge"],
                row["n_members"],
                row["neutral_mass"],
                row["mono_mz"],
                row["rt_width"],
                row["score"],
                row["quality"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def suspicious_low_score(conn, limit, max_score):
    print_header(f"Suspicious distributions: score <= {max_score:g}")

    rows = conn.execute(
        """
        SELECT
            distribution_id,
            charge,
            n_members,
            neutral_mass,
            mono_mz,
            rt_end - rt_start AS rt_width,
            score,
            quality,
            ms1_start,
            ms1_end
        FROM distributions
        WHERE score <= ?
        ORDER BY score ASC, quality ASC
        LIMIT ?
        """,
        (max_score, limit),
    ).fetchall()

    print_table(
        [
            "distribution_id",
            "z",
            "members",
            "mass",
            "mono_mz",
            "rt_width",
            "score",
            "quality",
            "ms1",
        ],
        [
            (
                row["distribution_id"],
                row["charge"],
                row["n_members"],
                row["neutral_mass"],
                row["mono_mz"],
                row["rt_width"],
                row["score"],
                row["quality"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def suspicious_high_charge_low_mass(conn, limit, max_mass, min_charge):
    print_header(f"Suspicious distributions: charge >= {min_charge} and mass <= {max_mass:g} Da")

    rows = conn.execute(
        """
        SELECT
            distribution_id,
            charge,
            n_members,
            neutral_mass,
            mono_mz,
            rt_end - rt_start AS rt_width,
            score,
            quality,
            ms1_start,
            ms1_end
        FROM distributions
        WHERE charge >= ?
          AND neutral_mass <= ?
        ORDER BY charge DESC, neutral_mass ASC
        LIMIT ?
        """,
        (min_charge, max_mass, limit),
    ).fetchall()

    print_table(
        [
            "distribution_id",
            "z",
            "members",
            "mass",
            "mono_mz",
            "rt_width",
            "score",
            "quality",
            "ms1",
        ],
        [
            (
                row["distribution_id"],
                row["charge"],
                row["n_members"],
                row["neutral_mass"],
                row["mono_mz"],
                row["rt_width"],
                row["score"],
                row["quality"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def suspicious_large_analytes(conn, limit, min_distributions):
    if not table_exists(conn, "analytes"):
        return

    print_header(f"Suspicious analytes: n_distributions >= {min_distributions}")

    rows = conn.execute(
        """
        SELECT
            analyte_id,
            neutral_mass,
            rt_start,
            rt_apex,
            rt_end,
            rt_end - rt_start AS rt_width,
            charge_min,
            charge_max,
            n_distributions,
            score,
            ms1_start,
            ms1_end
        FROM analytes
        WHERE n_distributions >= ?
        ORDER BY n_distributions DESC, score ASC
        LIMIT ?
        """,
        (min_distributions, limit),
    ).fetchall()

    print_table(
        [
            "analyte_id",
            "mass",
            "rt_width",
            "charges",
            "n_dists",
            "score",
            "ms1",
        ],
        [
            (
                row["analyte_id"],
                row["neutral_mass"],
                row["rt_width"],
                f"{row['charge_min']}-{row['charge_max']}",
                row["n_distributions"],
                row["score"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def most_complex_distributions(conn, limit):
    print_header("Most complex distributions by isotope member count")

    rows = conn.execute(
        """
        SELECT
            distribution_id,
            charge,
            n_members,
            neutral_mass,
            mono_mz,
            rt_end - rt_start AS rt_width,
            score,
            quality,
            ms1_start,
            ms1_end
        FROM distributions
        ORDER BY n_members DESC, quality DESC
        LIMIT ?
        """,
        (limit,),
    ).fetchall()

    print_table(
        [
            "distribution_id",
            "z",
            "members",
            "mass",
            "mono_mz",
            "rt_width",
            "score",
            "quality",
            "ms1",
        ],
        [
            (
                row["distribution_id"],
                row["charge"],
                row["n_members"],
                row["neutral_mass"],
                row["mono_mz"],
                row["rt_width"],
                row["score"],
                row["quality"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def highest_quality_distributions(conn, limit):
    print_header("Highest-quality distributions")

    rows = conn.execute(
        """
        SELECT
            distribution_id,
            charge,
            n_members,
            neutral_mass,
            mono_mz,
            rt_end - rt_start AS rt_width,
            score,
            quality,
            ms1_start,
            ms1_end
        FROM distributions
        ORDER BY quality DESC, score DESC
        LIMIT ?
        """,
        (limit,),
    ).fetchall()

    print_table(
        [
            "distribution_id",
            "z",
            "members",
            "mass",
            "mono_mz",
            "rt_width",
            "score",
            "quality",
            "ms1",
        ],
        [
            (
                row["distribution_id"],
                row["charge"],
                row["n_members"],
                row["neutral_mass"],
                row["mono_mz"],
                row["rt_width"],
                row["score"],
                row["quality"],
                f"{row['ms1_start']}-{row['ms1_end']}",
            )
            for row in rows
        ],
    )


def member_consistency(conn):
    if not table_exists(conn, "distribution_members"):
        return

    print_header("Distribution member consistency")

    total_distributions = count_rows(conn, "distributions") or 0

    mismatched = conn.execute(
        """
        SELECT count(*) AS n
        FROM (
            SELECT
                d.distribution_id,
                d.n_members AS expected,
                count(dm.feature_id) AS observed
            FROM distributions d
            LEFT JOIN distribution_members dm
              ON dm.distribution_id = d.distribution_id
            GROUP BY d.distribution_id
            HAVING expected != observed
        )
        """
    ).fetchone()["n"]

    duplicate_members = conn.execute(
        """
        SELECT count(*) AS n
        FROM (
            SELECT distribution_id, feature_id, count(*) AS n
            FROM distribution_members
            GROUP BY distribution_id, feature_id
            HAVING n > 1
        )
        """
    ).fetchone()["n"]

    orphan_members = conn.execute(
        """
        SELECT count(*) AS n
        FROM distribution_members dm
        LEFT JOIN distributions d
          ON d.distribution_id = dm.distribution_id
        WHERE d.distribution_id IS NULL
        """
    ).fetchone()["n"]

    rows = [
        ("total_distributions", total_distributions),
        ("mismatched_member_counts", mismatched),
        ("duplicate_member_rows", duplicate_members),
        ("orphan_member_rows", orphan_members),
    ]

    print_key_values(rows)


def analyte_consistency(conn):
    if not table_exists(conn, "analytes") or not table_exists(conn, "analyte_members"):
        return

    print_header("Analyte member consistency")

    total_analytes = count_rows(conn, "analytes") or 0

    mismatched = conn.execute(
        """
        SELECT count(*) AS n
        FROM (
            SELECT
                a.analyte_id,
                a.n_distributions AS expected,
                count(am.distribution_id) AS observed
            FROM analytes a
            LEFT JOIN analyte_members am
              ON am.analyte_id = a.analyte_id
            GROUP BY a.analyte_id
            HAVING expected != observed
        )
        """
    ).fetchone()["n"]

    orphan_members = conn.execute(
        """
        SELECT count(*) AS n
        FROM analyte_members am
        LEFT JOIN analytes a
          ON a.analyte_id = am.analyte_id
        WHERE a.analyte_id IS NULL
        """
    ).fetchone()["n"]

    rows = [
        ("total_analytes", total_analytes),
        ("mismatched_member_counts", mismatched),
        ("orphan_member_rows", orphan_members),
    ]

    print_key_values(rows)


def write_json_report(path, report):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def json_summary(conn, db_path, params):
    summary = {
        "database": str(db_path),
        "size_bytes": Path(db_path).stat().st_size,
        "parameters": params,
        "tables": {},
        "quantiles": {},
    }

    for table, n in table_counts(conn):
        summary["tables"][table] = n

    if table_exists(conn, "distributions"):
        summary["quantiles"]["distribution_score"] = quantiles(
            fetch_values(conn, "SELECT score FROM distributions")
        )
        summary["quantiles"]["distribution_quality"] = quantiles(
            fetch_values(conn, "SELECT quality FROM distributions")
        )
        summary["quantiles"]["distribution_rt_width"] = quantiles(
            fetch_values(conn, "SELECT rt_end - rt_start FROM distributions")
        )
        summary["quantiles"]["distribution_neutral_mass"] = quantiles(
            fetch_values(conn, "SELECT neutral_mass FROM distributions")
        )

    if table_exists(conn, "features"):
        summary["quantiles"]["feature_n_points"] = quantiles(
            fetch_values(conn, "SELECT n_points FROM features")
        )
        summary["quantiles"]["feature_rt_width"] = quantiles(
            fetch_values(conn, "SELECT rt_end - rt_start FROM features")
        )

    if table_exists(conn, "analytes"):
        summary["quantiles"]["analyte_n_distributions"] = quantiles(
            fetch_values(conn, "SELECT n_distributions FROM analytes")
        )

    return summary


def inspect(args):
    conn = connect(args.sqlite)

    try:
        params = fetch_parameters(conn)

        print_header("Database metadata")
        print_key_values(metadata_rows(args.sqlite, params))

        print_header("Table counts")
        print_table(["table", "count"], table_counts(conn))

        if not table_exists(conn, "distributions"):
            raise SystemExit("missing table: distributions")

        member_consistency(conn)
        analyte_consistency(conn)

        print_distribution_by_charge(conn)
        print_distribution_by_members(conn)

        if args.charge_members:
            print_distribution_charge_members(conn)

        print_distribution_score_summaries(conn)

        if args.features:
            print_feature_summaries(conn)

        if args.analytes:
            print_analyte_summaries(conn)

        highest_quality_distributions(conn, args.limit)
        most_complex_distributions(conn, args.limit)

        suspicious_tiny_rt(conn, args.limit, args.tiny_rt)
        suspicious_huge_rt(conn, args.limit, args.huge_rt)
        suspicious_low_score(conn, args.limit, args.low_score)
        suspicious_high_charge_low_mass(conn, args.limit, args.low_mass, args.high_charge)
        suspicious_large_analytes(conn, args.limit, args.large_analyte)

        if args.json_out:
            report = json_summary(conn, args.sqlite, params)
            write_json_report(args.json_out, report)
            print()
            print(f"wrote JSON summary: {args.json_out}")

    finally:
        conn.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Inspect a distribution deconvolution SQLite file."
    )

    parser.add_argument("sqlite", type=Path)

    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--json-out", type=Path)

    parser.add_argument("--tiny-rt", type=float, default=0.02)
    parser.add_argument("--huge-rt", type=float, default=2.0)
    parser.add_argument("--low-score", type=float, default=0.35)
    parser.add_argument("--low-mass", type=float, default=500.0)
    parser.add_argument("--high-charge", type=int, default=4)
    parser.add_argument("--large-analyte", type=int, default=10)

    parser.add_argument("--features", action="store_true")
    parser.add_argument("--analytes", action="store_true")
    parser.add_argument("--charge-members", action="store_true")

    return parser.parse_args()


def main():
    args = parse_args()
    inspect(args)


if __name__ == "__main__":
    main()
