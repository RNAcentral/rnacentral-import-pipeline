import groovy.transform.Canonical

@Canonical
class DataSource {
  String name;
  List<Input> inputs;
  Process process;

  static DataSource build(db_name, database) {
    def index = 0
    def inputs = []
    database.inputs.each { input_entry ->
      inputs << Input.build(db_name, index, input_entry.key, input_entry.value)
      index++
    }

    def group = database.get('group', 'external')
    return new DataSource([
      name: db_name,
      inputs: inputs,
      process: [
        command: database.get('command', "rnac $group $db_name"),
        produces: database.get('produces', '*.csv'),
        directives: [
          memory: database.get('memory', 2.GB),
          tag: db_name,
          group: "/rnc/process/$db_name"
        ]
      ]
    ]);
  }

  def script(input_files) {
    String prefix = ''
    List<String> arguments = []
    int gzip_count = inputs.inject(0) { a, i -> a + (i.produces.endsWith('.gz') ? 1 : 0) }
    if (gzip_count == 0) {
      arguments.addAll(input_files)
    } else if (gzip_count == 1) {
      prefix = "zcat *.gz | "
      arguments = input_files.inject([]) { acc, fn -> acc << (fn.endsWith('.gz') ? '-' : fn) }
    } else {
      prefix = "gzip -d *.gz\n"
      arguments = input_files.inject([]) { acc, fn -> fn.replace('.gz', '') }
    }
    return "$prefix $process.command ${arguments.join(' ')}".trim();
  }

  Boolean is_parallel_task() {
    return this.inputs.any { input -> input.produces.contains('*') }
  }
}

@Canonical
class Input {
  String source;
  int index;
  String command;
  List<String> arguments;
  List<String> exclude;
  String produces;
  Map directives;

  static Map defaults(db_name, index, input_name) {
    return [
      source: db_name,
      index: index,
      command: 'fetch generic',
      exclude: [],
      directives: [
        memory: 2.GB,
        tag: "$db_name $input_name",
        group: "/rnac/fetch/$db_name/$input_name"
      ],
    ]
  }

  static Input from_string(db_name, index, input_name, input) {
    def f = new File(input)
    def produces = f.getName()
    def spec = [arguments: [input], exclude: [], produces: produces];

    return new Input(Input.defaults(db_name, index, input_name) + spec);
  }

  static Input from_mysql(db_name, index, input_name, mysql) {
    def options = mysql
      .subMap(['host', 'port', 'user', 'database'])
      .inject([]) { result, entry -> result << "--$entry.key" << entry.value }
      .join(' ')

    def spec = [
      command: "mysql ${options}",
      arguments: [options.query],
      produces: mysql.produces,
    ];

    return new Input(Input.defaults(db_name, index, input_name) + spec);
  }

  static Input url_file(db_name, index, input_name, input) {

      def url_file = file("$db_name-urls-${index}.txt")
      url_file.text = input.url_file.join('\n')
      def spec = [
        command: input.get('command', 'fetch url-file'),
        arguments: [url_file],
        produces: input.produces,
      ];

      return new Input(Input.defaults(db_name, index, input_name) + spec);
  }

  static Input build(db_name, index, input_name, input) {

    if (input instanceof String) {
      return Input.from_string(db_name, index, input_name, input);
    }

    if (input.containsKey('command') && input.command == 'mysql') {
      return Input.from_mysql(input);
    }

    if (input.containsKey('url_file')) {
      return Input.url_file(db_name, index, input_name, input);
    }

    Map spec = Input.defaults(db_name, index, input_name);
    spec.arguments = input.remote instanceof List ? input.remote : [input.remote];
    spec.command = input.get('command', spec.command);
    spec.directives.memory = input.get('memory', spec.directives.memory);

    String produces;
    if (input.containsKey('produces')) {
      spec.produces = input.produces
    } else {
      def f = new File(input.remote)
      spec.produces = f.getName()
    }

    return new Input(spec);
  }

  def script() {
    def args = arguments.inject([]) { agg, a -> agg << "'$a'" }
    if (command.startsWith('mysql')) {
      return "$command < ${args[0]} > '$produces'"
    }
    return "$command ${args.join(' ')} '$produces'"
  }
}

@Canonical
class Process {
  String command;
  String produces;
  Map directives;
}
