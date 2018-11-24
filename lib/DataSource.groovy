class DataSource {
  static Map build(db_name, database) {
    def inputs = database.get('inputs', []).inject([]) { acc, input_entry ->
      acc << Input.build(db_name, acc.length + 1, input_entry.key, input_entry.value)
    }

    Map process = Parse.build(db_name, database.process);
    return [
      name: db_name,
      inputs: inputs,
      process: process
    ];
  }

  static def process_script(Map source, List<String> input_files) {
    String prefix = ''
    List<String> arguments = []
    int gzip_count = source.inputs.inject(0) { a, i -> a + (i.produces.endsWith('.gz') ? 1 : 0) }
    if (gzip_count == 0) {
      arguments.addAll(input_files)
    } else if (gzip_count == 1 && !source.process.force_decompress) {
      prefix = "zcat *.gz | "
      arguments = input_files.inject([]) { acc, fn -> acc << (fn.endsWith('.gz') ? '-' : fn) }
    } else {
      prefix = "gzip -fd *.gz\n"
      arguments = input_files.inject([]) { acc, fn -> acc << fn.replace('.gz', '') }
    }
    return "$prefix $source.process.command ${arguments.join(' ')}".trim();
  }

  static Boolean is_parallel_task(Map source) {
    return source.inputs.any { input -> input.produces.contains('*') }
  }
}
