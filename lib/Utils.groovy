static def orderedScripts(names) {
  return names.sort(false) { raw ->
    def parts = raw.split("__");
    return parts[0].toInteger();
  };
}

static def write_ordered(output, scriptNames) {
  output.withWriter { t ->
    orderedScripts(scriptNames).each { s ->
      t << s << "\n"
    }
  }

  return output;
}

static def must_release(to_import, databases) {
  return databases.inject(false) { s, e ->
    s || (e.value && databases[e.key].get('release', true))
  };
}

// Flags derived from params.databases. Shared because main.nf and
// import-data.nf are both entry points and need the same values.
static def ENSEMBL_DIVISIONS() {
  return ['vertebrates', 'plants', 'fungi', 'protists', 'metazoa'];
}

// Ensembl is special-cased: it runs if any of its divisions do.
static def will_run(key, db) {
  if (!(db instanceof Map)) {
    return false;
  }
  if (key == 'ensembl') {
    return ENSEMBL_DIVISIONS().any { d -> db[d]?.get('run', false) };
  }
  return db.get('run', false);
}

static def ensembl_runs(databases) {
  def ensembl = databases?.ensembl;
  return ensembl != null && ENSEMBL_DIVISIONS().any { d -> ensembl[d]?.get('run', false) };
}

static def should_release(databases) {
  return databases.any { key, db -> will_run(key, db) && db.get('release', true) };
}

static def needs_publications(databases, skip) {
  return !skip && databases.any { key, db -> will_run(key, db) };
}

static def needs_taxonomy(databases) {
  return databases.any { key, db -> will_run(key, db) && db.get('needs_taxonomy', false) };
}
