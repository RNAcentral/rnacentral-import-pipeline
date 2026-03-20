# -*- coding: utf-8 -*-


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
