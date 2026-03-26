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

    conn = psycopg2.connect(conn_str)
    registered = []
    try:
        with conn.cursor() as cur:
            for job_id in ids:
                cur.execute(
                    """
                    INSERT INTO litscan_job (job_id, display_id, submitted, status)
                    VALUES (%s, %s, %s, 'pending')
                    ON CONFLICT (job_id) DO NOTHING
                    """,
                    (job_id.lower(), job_id, datetime.datetime.now()),
                )
                if cur.rowcount > 0:
                    registered.append(job_id)
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
