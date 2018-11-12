class DataSource {
  static Map defaults(db_name) {
    return [
      name: db_name,
      inputs: [],
      process: [
        command: "rnac external $db_name",
        produces: '*.csv',
        directives: [
          memory: 2.GB,
          tag: db_name,
          group: "/rnc/process/$db_name"
        ]
      ]
    ];
  }

  static Map build(db_name, database) {
    def index = 0
    def inputs = []
    database.inputs.each { input_entry ->
      inputs << Input.build(db_name, index, input_entry.key, input_entry.value)
      index++
    }

    Map update = database.subMap('command', 'produces', 'memory');
    update.inputs = inputs;
    return DataSource.defaults(db_name) + update;
  }

  static def process_script(Map source, List<String> input_files) {
    String prefix = ''
    List<String> arguments = []
    int gzip_count = source.inputs.inject(0) { a, i -> a + (i.produces.endsWith('.gz') ? 1 : 0) }
    if (gzip_count == 0) {
      arguments.addAll(input_files)
    } else if (gzip_count == 1) {
      prefix = "zcat *.gz | "
      arguments = input_files.inject([]) { acc, fn -> acc << (fn.endsWith('.gz') ? '-' : fn) }
    } else {
      prefix = "gzip -d *.gz\n"
      arguments = input_files.inject([]) { acc, fn -> fn.replace('.gz', '') }
    }
    return "$prefix $source.process.command ${arguments.join(' ')}".trim();
  }

  static Boolean is_parallel_task(Map source) {
    return source.inputs.any { input -> input.produces.contains('*') }
  }
}
