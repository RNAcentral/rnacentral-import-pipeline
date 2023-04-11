# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import pickle
import threading


def pickle_stream(stream, handle, *args, **kwargs):
    for entry in stream:
        pickle.dump(entry, handle, *args, **kwargs)


def unpickle_stream(handle, *args, **kwargs):
    try:
        while True:
            yield pickle.load(handle)
    except EOFError:
        return


## From https://stackoverflow.com/a/46723144/3249000 - make the cache store async objects properly
class Cacheable:
    def __init__(self, co):
        self.co = co
        self.done = False
        self.result = None
        self.lock = threading.RLock()
        ## This needs to be a re-rntrant lock so it is only release by the coroutine that acquired it

    def __await__(self):
        with self.lock:
            if self.done:
                return self.result
            self.result = yield from self.co.__await__()
            self.done = True
            return self.result

def cacheable(f):
    def wrapped(*args, **kwargs):
        r = f(*args, **kwargs)
        return Cacheable(r)
    return wrapped

