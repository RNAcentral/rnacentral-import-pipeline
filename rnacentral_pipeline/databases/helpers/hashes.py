# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import hashlib

try:
    import crcmod
except ImportError:
    crcmod = None


def md5(data):
    """
    Get the MD5 hash as a string.
    """
    return hashlib.md5(data).hexdigest()


_POLY64REVh = 0xD8000000


def _build_crc64_tables():
    table_h = [0] * 256
    table_l = [0] * 256
    for i in range(256):
        partl = i
        parth = 0
        for _ in range(8):
            rflag = partl & 1
            partl >>= 1
            if parth & 1:
                partl |= 1 << 31
            parth >>= 1
            if rflag:
                parth ^= _POLY64REVh
        table_h[i] = parth
        table_l[i] = partl
    return table_h, table_l


_CRC64_TABLE_H, _CRC64_TABLE_L = _build_crc64_tables()


def crc64_python(input_string):
    """
    Python re-implementation of SWISS::CRC64
    Adapted from:
    http://code.activestate.com/recipes/259177-crc64-calculate-the-cyclic-redundancy-check/
    The generator polynomial is x64 + x4 + x3 + x + 1.

    This is kept as a fallback incase crcmod is not available. It is quite a bit slower though.

    Implementation below tested against 49,991 samples from ENA and gave identical results

    """
    crcl = 0
    crch = 0
    for item in input_string:
        shr = 0
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        tableindex = (crcl ^ ord(item)) & 0xFF

        crch = temp1h ^ _CRC64_TABLE_H[tableindex]
        crcl = temp1l ^ _CRC64_TABLE_L[tableindex]
    return "%08X%08X" % (crch, crcl)


if crcmod is not None:
    _crc64_swiss = crcmod.mkCrcFun(0x1000000000000001B, initCrc=0, rev=True, xorOut=0)
else:
    _crc64_swiss = None


def crc64(input_string):
    if _crc64_swiss is None:
        if isinstance(input_string, bytes):
            input_string = input_string.decode("ascii")
        return crc64_python(input_string)
    if isinstance(input_string, str):
        input_string = input_string.encode("ascii")
    return "%016X" % _crc64_swiss(input_string)
