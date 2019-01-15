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
