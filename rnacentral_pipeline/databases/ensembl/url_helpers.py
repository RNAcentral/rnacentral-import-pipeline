# -*- coding: utf-8 -*-

from fnmatch import fnmatch
from ftplib import FTP
from posixpath import basename, dirname
from urllib.parse import urlparse


def lowercase_assembly_in_url(url: str) -> str:
    prefix, separator, filename = url.rpartition("/")
    if not separator:
        filename = prefix
        prefix = ""

    parts = filename.split(".", 2)
    if len(parts) < 3:
        return url

    parts[1] = parts[1].lower()
    lowered = ".".join(parts)
    if not separator:
        return lowered
    return f"{prefix}/{lowered}"


def resolve_ftp_urls(url: str) -> list[str]:
    if "*" not in url:
        return [url]

    for candidate in [url, lowercase_assembly_in_url(url)]:
        matches = _resolve_ftp_glob(candidate)
        if matches:
            return matches

    return []


def _resolve_ftp_glob(url: str) -> list[str]:
    parsed = urlparse(url)
    if parsed.scheme != "ftp":
        return [url]

    remote_dir = dirname(parsed.path)
    pattern = basename(parsed.path)
    with FTP(parsed.hostname) as ftp:
        ftp.login()
        names = ftp.nlst(remote_dir)

    return [
        f"ftp://{parsed.hostname}{remote_dir}/{basename(name)}"
        for name in names
        if fnmatch(basename(name), pattern)
    ]
