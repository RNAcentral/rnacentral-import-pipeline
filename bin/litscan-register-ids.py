#!/usr/bin/env python3
"""
Register new RNA IDs in the litscan_job table for scanning.

Usage:
    litscan-register-ids.py <input_file> <output_file>

Reads IDs from <input_file> (one per line), inserts each into litscan_job
with status='pending' using ON CONFLICT DO NOTHING, and writes the
successfully registered IDs to <output_file>.

Environment variables:
    PSYCOPG_CONN   PostgreSQL connection URI
"""
import datetime
import os
import sys

import psycopg2
from psycopg2.extras import execute_values


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]
    conn_str = os.environ["PSYCOPG_CONN"]

    with open(input_file) as fh:
        ids = [line.strip() for line in fh if line.strip()]

    if not ids:
        open(output_file, "w").close()
        return

    now = datetime.datetime.now()
    rows = [(job_id.lower(), job_id, now) for job_id in ids]

    conn = psycopg2.connect(conn_str)
    registered = []
    try:
        with conn.cursor() as cur:
            execute_values(
                cur,
                """
                INSERT INTO litscan_job (job_id, display_id, submitted, status)
                SELECT v.job_id, v.display_id, v.submitted, 'pending'
                FROM (VALUES %s) AS v(job_id, display_id, submitted)
                WHERE NOT EXISTS (
                    SELECT 1 FROM litscan_job j
                    WHERE j.job_id = v.job_id
                    AND j.status IN ('pending', 'success')
                )
                RETURNING job_id
                """,
                rows,
                page_size=1000,
            )
            registered = [r[0] for r in cur.fetchall()]
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()

    with open(output_file, "w") as fh:
        for job_id in registered:
            fh.write(job_id + "\n")

    print(f"Registered {len(registered)} of {len(ids)} IDs")


if __name__ == "__main__":
    main()
