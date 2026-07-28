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
static def will_run(db) {
  return db instanceof Map && db.get('run', false);
}

static def should_release(databases) {
  return databases.any { _key, db -> will_run(db) && db.get('release', true) };
}

static def needs_publications(databases, skip) {
  return !skip && databases.any { _key, db -> will_run(db) };
}

static def needs_taxonomy(databases) {
  return databases.any { _key, db -> will_run(db) && db.get('needs_taxonomy', false) };
}
