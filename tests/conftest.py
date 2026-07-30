# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

"""
Root test configuration.

The unit test suite must run without network access: any test that needs data
from an external service is expected to mock it (see, for example,
``tests/databases/pdb/conftest.py``). To enforce that, this module blocks all
outbound socket connections to non-local addresses during the test run. A test
that genuinely needs live network access must be marked ``@pytest.mark.network``;
those tests are skipped unless network access is explicitly enabled by setting
the ``RNAC_TEST_ALLOW_NETWORK`` environment variable.

This guard means a newly added (or regressed) unmocked network call fails loudly
with ``NetworkAccessBlocked`` instead of silently depending on the outside world.
"""

import json
import os
import socket
from pathlib import Path

import pytest

from tests import _cassette

ALLOW_NETWORK = bool(os.environ.get("RNAC_TEST_ALLOW_NETWORK"))

# Hosts we always allow - loopback is needed for things like local databases.
_ALLOWED_HOSTS = {"127.0.0.1", "::1", "localhost", "0.0.0.0"}

_real_socket_connect = socket.socket.connect
_real_create_connection = socket.create_connection
_real_getaddrinfo = socket.getaddrinfo


class NetworkAccessBlocked(RuntimeError):
    """Raised when a test attempts a non-local network connection."""


def _host_of(address):
    if isinstance(address, (tuple, list)) and address:
        return address[0]
    return address


def _is_allowed(host):
    return host in _ALLOWED_HOSTS


def _blocked(host):
    raise NetworkAccessBlocked(
        f"Blocked network access to {host!r} during tests. Mock the external "
        "call (see tests/databases/pdb/conftest.py for the pattern) or mark the "
        "test with @pytest.mark.network if it genuinely needs live access."
    )


def _guarded_connect(self, address, *args, **kwargs):
    host = _host_of(address)
    if not _is_allowed(host):
        _blocked(host)
    return _real_socket_connect(self, address, *args, **kwargs)


def _guarded_create_connection(address, *args, **kwargs):
    host = _host_of(address)
    if not _is_allowed(host):
        _blocked(host)
    return _real_create_connection(address, *args, **kwargs)


def _guarded_getaddrinfo(host, *args, **kwargs):
    if not _is_allowed(host):
        _blocked(host)
    return _real_getaddrinfo(host, *args, **kwargs)


def pytest_configure(config):
    _cassette.install(ALLOW_NETWORK)
    _install_phylogeny_cache()
    if ALLOW_NETWORK:
        return
    socket.socket.connect = _guarded_connect
    socket.create_connection = _guarded_create_connection
    socket.getaddrinfo = _guarded_getaddrinfo


def pytest_unconfigure(config):
    _cassette.uninstall()
    _restore_phylogeny_cache()
    socket.socket.connect = _real_socket_connect
    socket.create_connection = _real_create_connection
    socket.getaddrinfo = _real_getaddrinfo


def pytest_runtest_setup(item):
    """Skip tests marked ``network`` unless live access is explicitly enabled."""
    if item.get_closest_marker("network") and not ALLOW_NETWORK:
        pytest.skip("requires network access (set RNAC_TEST_ALLOW_NETWORK to run)")


# ---------------------------------------------------------------------------
# Cached ENA/UniProt taxonomy lookups
# ---------------------------------------------------------------------------
# A large number of database parser tests resolve a taxon id to its lineage /
# scientific name through rnacentral_pipeline.databases.helpers.phylogeny, which
# hits the EBI ENA taxonomy REST API. Rather than mock each test individually we
# back the shared phylogeny helpers with a vendored, on-disk cache of real API
# responses so the tests run offline.
#
# To refresh or extend the cache, run the suite once with live access:
#   RNAC_TEST_ALLOW_NETWORK=1 uv run pytest <tests that need new taxids>
# any taxon id / species not already cached is fetched from the real API and
# written back into tests/data/phylogeny/.

_PHYLO_DIR = Path(__file__).parent / "data" / "phylogeny"
_PHYLO_TAX_CACHE = _PHYLO_DIR / "ena_taxonomy.json"
_PHYLO_NAME_CACHE = _PHYLO_DIR / "ena_taxid.json"


def _load_cache(path):
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _store_cache(path, data):
    _PHYLO_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True))


# Modules that import ``taxid`` by name and so need their own binding patched.
_TAXID_REBOUND_MODULES = (
    "rnacentral_pipeline.rnacentral.r2dt.models.ribovision",
    "rnacentral_pipeline.rnacentral.r2dt.models.crw",
)

# Saved original callables so they can be restored in pytest_unconfigure.
_phylo_originals = {}


def _install_phylogeny_cache():
    """
    Replace the network-backed phylogeny helpers with cache-backed versions.

    This is installed in ``pytest_configure`` rather than a fixture so that the
    patch is in place before any module/session scoped fixture (which may resolve
    taxonomy during its own setup) runs.
    """
    from rnacentral_pipeline.databases.helpers import phylogeny as phy

    tax_cache = _load_cache(_PHYLO_TAX_CACHE)
    name_cache = _load_cache(_PHYLO_NAME_CACHE)

    real_phylogeny = getattr(phy.phylogeny, "__wrapped__", phy.phylogeny)
    real_taxid = getattr(phy.taxid, "__wrapped__", phy.taxid)

    def fake_phylogeny(taxon_id):
        key = str(taxon_id)
        if key in tax_cache:
            return tax_cache[key]
        if ALLOW_NETWORK:
            data = real_phylogeny(taxon_id)
            tax_cache[key] = data
            _store_cache(_PHYLO_TAX_CACHE, tax_cache)
            return data
        raise phy.FailedTaxonId(
            f"No cached phylogeny for taxon id {taxon_id}. Record it with "
            "RNAC_TEST_ALLOW_NETWORK=1."
        )

    def fake_taxid(species):
        key = str(species)
        if key in name_cache:
            return name_cache[key]
        if ALLOW_NETWORK:
            value = real_taxid(species)
            name_cache[key] = value
            _store_cache(_PHYLO_NAME_CACHE, name_cache)
            return value
        raise phy.FailedTaxonId(
            f"No cached taxid for species {species!r}. Record it with "
            "RNAC_TEST_ALLOW_NETWORK=1."
        )

    _phylo_originals[(phy, "phylogeny")] = phy.phylogeny
    _phylo_originals[(phy, "taxid")] = phy.taxid
    phy.phylogeny = fake_phylogeny
    phy.taxid = fake_taxid

    for module_path in _TAXID_REBOUND_MODULES:
        try:
            module = __import__(module_path, fromlist=["taxid"])
        except Exception:
            continue
        if hasattr(module, "taxid"):
            _phylo_originals[(module, "taxid")] = module.taxid
            module.taxid = fake_taxid


def _restore_phylogeny_cache():
    for (obj, name), original in _phylo_originals.items():
        setattr(obj, name, original)
    _phylo_originals.clear()
