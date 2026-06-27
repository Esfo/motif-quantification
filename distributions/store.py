import json
import sqlite3
from pathlib import Path


def connect_db(path, overwrite=False):
    path = Path(path)

    if overwrite and path.exists():
        path.unlink()

    path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA foreign_keys=ON")

    return conn


def init_schema(conn, store_edges=False):
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS parameters (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS scans (
            ms1_index INTEGER PRIMARY KEY,
            spectrum_index INTEGER NOT NULL,
            scan_id TEXT NOT NULL,
            rt REAL NOT NULL,
            tic REAL NOT NULL,
            n_points INTEGER NOT NULL
        );

        CREATE TABLE IF NOT EXISTS lines (
            line_id INTEGER PRIMARY KEY,
            mz_mean REAL NOT NULL,
            mz_min REAL NOT NULL,
            mz_max REAL NOT NULL,
            rt_start REAL NOT NULL,
            rt_end REAL NOT NULL,
            ms1_start INTEGER NOT NULL,
            ms1_end INTEGER NOT NULL,
            n_points INTEGER NOT NULL
        );

        CREATE TABLE IF NOT EXISTS features (
            feature_id INTEGER PRIMARY KEY,
            line_id INTEGER NOT NULL,
            mz_mean REAL NOT NULL,
            mz_min REAL NOT NULL,
            mz_max REAL NOT NULL,
            rt_start REAL NOT NULL,
            rt_apex REAL NOT NULL,
            rt_end REAL NOT NULL,
            ms1_start INTEGER NOT NULL,
            ms1_apex INTEGER NOT NULL,
            ms1_end INTEGER NOT NULL,
            height REAL NOT NULL,
            area REAL NOT NULL,
            n_points INTEGER NOT NULL,
            quality REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS feature_traces (
            feature_id INTEGER PRIMARY KEY,
            n INTEGER NOT NULL,
            scans BLOB NOT NULL,
            rts BLOB NOT NULL,
            mzs BLOB NOT NULL,
            intensities BLOB NOT NULL
        );

        CREATE TABLE IF NOT EXISTS distributions (
            distribution_id INTEGER PRIMARY KEY,
            charge INTEGER NOT NULL,
            neutral_mass REAL NOT NULL,
            mono_mz REAL NOT NULL,
            rt_start REAL NOT NULL,
            rt_apex REAL NOT NULL,
            rt_end REAL NOT NULL,
            ms1_start INTEGER NOT NULL,
            ms1_apex INTEGER NOT NULL,
            ms1_end INTEGER NOT NULL,
            n_members INTEGER NOT NULL,
            score REAL NOT NULL,
            quality REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS distribution_members (
            distribution_id INTEGER NOT NULL,
            feature_id INTEGER NOT NULL,
            isotope_index INTEGER NOT NULL,
            member_score REAL NOT NULL,
            PRIMARY KEY (distribution_id, feature_id)
        );

        CREATE TABLE IF NOT EXISTS analytes (
            analyte_id INTEGER PRIMARY KEY,
            neutral_mass REAL NOT NULL,
            rt_start REAL NOT NULL,
            rt_apex REAL NOT NULL,
            rt_end REAL NOT NULL,
            ms1_start INTEGER NOT NULL,
            ms1_apex INTEGER NOT NULL,
            ms1_end INTEGER NOT NULL,
            charge_min INTEGER NOT NULL,
            charge_max INTEGER NOT NULL,
            n_distributions INTEGER NOT NULL,
            score REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS analyte_members (
            analyte_id INTEGER NOT NULL,
            distribution_id INTEGER NOT NULL,
            charge INTEGER NOT NULL,
            PRIMARY KEY (analyte_id, distribution_id)
        );
        """
    )

    if store_edges:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS isotope_edges (
                edge_id INTEGER PRIMARY KEY,
                left_feature_id INTEGER NOT NULL,
                right_feature_id INTEGER NOT NULL,
                charge INTEGER NOT NULL,
                isotope_step INTEGER NOT NULL,
                mz_error REAL NOT NULL,
                mz_error_ppm REAL NOT NULL,
                rt_shift REAL NOT NULL,
                rt_overlap REAL NOT NULL,
                intensity_score REAL NOT NULL,
                score REAL NOT NULL
            );
            """
        )


def create_indexes(conn, store_edges=False):
    conn.executescript(
        """
        CREATE INDEX IF NOT EXISTS idx_scans_rt ON scans(rt);

        CREATE INDEX IF NOT EXISTS idx_lines_mz ON lines(mz_mean);
        CREATE INDEX IF NOT EXISTS idx_lines_rt ON lines(rt_start, rt_end);
        CREATE INDEX IF NOT EXISTS idx_lines_ms1 ON lines(ms1_start, ms1_end);

        CREATE INDEX IF NOT EXISTS idx_features_mz ON features(mz_mean);
        CREATE INDEX IF NOT EXISTS idx_features_rt ON features(rt_apex);
        CREATE INDEX IF NOT EXISTS idx_features_window ON features(rt_start, rt_end);
        CREATE INDEX IF NOT EXISTS idx_features_ms1 ON features(ms1_start, ms1_end);

        CREATE INDEX IF NOT EXISTS idx_distributions_mass ON distributions(neutral_mass);
        CREATE INDEX IF NOT EXISTS idx_distributions_rt ON distributions(rt_apex);
        CREATE INDEX IF NOT EXISTS idx_distributions_window ON distributions(rt_start, rt_end);
        CREATE INDEX IF NOT EXISTS idx_distributions_ms1 ON distributions(ms1_start, ms1_end);
        CREATE INDEX IF NOT EXISTS idx_distributions_charge ON distributions(charge);

        CREATE INDEX IF NOT EXISTS idx_distribution_members_dist ON distribution_members(distribution_id);
        CREATE INDEX IF NOT EXISTS idx_distribution_members_feature ON distribution_members(feature_id);

        CREATE INDEX IF NOT EXISTS idx_analytes_mass ON analytes(neutral_mass);
        CREATE INDEX IF NOT EXISTS idx_analytes_rt ON analytes(rt_apex);
        CREATE INDEX IF NOT EXISTS idx_analytes_window ON analytes(rt_start, rt_end);
        CREATE INDEX IF NOT EXISTS idx_analyte_members_analyte ON analyte_members(analyte_id);
        CREATE INDEX IF NOT EXISTS idx_analyte_members_distribution ON analyte_members(distribution_id);
        """
    )

    if store_edges:
        conn.executescript(
            """
            CREATE INDEX IF NOT EXISTS idx_edges_left ON isotope_edges(left_feature_id);
            CREATE INDEX IF NOT EXISTS idx_edges_right ON isotope_edges(right_feature_id);
            CREATE INDEX IF NOT EXISTS idx_edges_charge ON isotope_edges(charge);
            """
        )


def write_parameters(conn, params):
    rows = [(key, json.dumps(value)) for key, value in sorted(params.items())]

    conn.executemany(
        "INSERT OR REPLACE INTO parameters(key, value) VALUES (?, ?)",
        rows,
    )


def write_rows(conn, table, rows, batch_size=50000):
    rows = iter(rows)

    first = None

    for first in rows:
        break

    if first is None:
        return

    fields = list(first)
    placeholders = ",".join("?" for _ in fields)
    field_sql = ",".join(fields)
    sql = f"INSERT INTO {table} ({field_sql}) VALUES ({placeholders})"

    batch = [[first[field] for field in fields]]

    for row in rows:
        batch.append([row[field] for field in fields])

        if len(batch) >= batch_size:
            conn.executemany(sql, batch)
            batch.clear()

    if batch:
        conn.executemany(sql, batch)
