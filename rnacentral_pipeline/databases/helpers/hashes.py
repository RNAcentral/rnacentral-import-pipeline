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


def md5(data):
    """
    Get the MD5 hash as a string.
    """
    return hashlib.md5(data).hexdigest()


def crc64(input_string):
    """
    Python re-implementation of SWISS::CRC64
    Adapted from:
    http://code.activestate.com/recipes/259177-crc64-calculate-the-cyclic-redundancy-check/
    """

    POLY64REVh = 0xd8000000
    CRCTableh = [0] * 256
    CRCTablel = [0] * 256
    isInitialized = False
    crcl = 0
    crch = 0
    if isInitialized is not True:
        isInitialized = True
        for i in range(256):
            partl = i
            parth = 0
            for _ in range(8):
                rflag = partl & 1
                partl >>= 1
                if parth & 1:
                    partl |= (1 << 31)
                parth >>= 1
                if rflag:
                    parth ^= POLY64REVh
            CRCTableh[i] = parth
            CRCTablel[i] = partl

    for item in input_string:
        shr = 0
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        tableindex = (crcl ^ ord(item)) & 0xFF

        crch = temp1h ^ CRCTableh[tableindex]
        crcl = temp1l ^ CRCTablel[tableindex]
    return "%08X%08X" % (crch, crcl)
