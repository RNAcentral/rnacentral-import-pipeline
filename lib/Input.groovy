import static nextflow.Nextflow.file

class Input {
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

  static Map from_string(db_name, index, input_name, input) {
    def f = new File(input)
    def produces = f.getName()
    def spec = [arguments: [input], produces: produces];

    return Input.defaults(db_name, index, input_name) + spec;
  }

  static Map from_mysql(db_name, index, input_name, mysql) {
    def options = mysql
      .subMap(['host', 'port', 'user', 'database'])
      .inject([]) { result, entry -> result << "--$entry.key" << entry.value }
      .join(' ')

    def query = new File(mysql.query);
    def produces = query.getBaseName() + '.tsv'
    if (mysql.containsKey('produces')) {
      produces = mysql.produces
    }

    Map spec = [
      command: "mysql ${options}",
      arguments: [file(mysql.query)],
      produces: produces,
    ];

    return Input.defaults(db_name, index, input_name) + spec;
  }

  static Map urls(db_name, index, input_name, input) {

      def url_file = file("$db_name-urls-${index}.txt")
      url_file.text = input.urls.join('\n')
      def spec = [
        command: input.get('command', 'fetch url-file'),
        arguments: [url_file],
        produces: input.produces,
      ];

      return Input.defaults(db_name, index, input_name) + spec;
  }

  static Map url_file(db_name, index, input_name, input) {

      def spec = [
        command: input.get('command', 'fetch url-file'),
        arguments: [file(input.url_file)],
        produces: input.produces,
      ];

      return Input.defaults(db_name, index, input_name) + spec;
  }

  static Map build(db_name, index, input_name, input) {

    if (input instanceof String) {
      return Input.from_string(db_name, index, input_name, input);
    }

    if (input.containsKey('command') && input.command == 'mysql') {
      return Input.from_mysql(db_name, index, input_name, input);
    }

    if (input.containsKey('urls')) {
      return Input.urls(db_name, index, input_name, input);
    }

    if (input.containsKey('url_file')) {
      return Input.url_file(db_name, index, input_name, input);
    }

    Map spec = Input.defaults(db_name, index, input_name);
    Map update = input.subMap('command', 'memory', 'exclude');
    update.arguments = input.remote instanceof List ? input.remote : [input.remote];
    spec += update

    String produces;
    if (input.containsKey('produces')) {
      spec.produces = input.produces
    } else {
      def f = new File(input.remote)
      spec.produces = f.getName()
    }

    return spec;
  }

  static String script(Map input) {
    def args = input.arguments.inject([]) { agg, a -> agg << "'$a'" }
    if (input.command.startsWith('mysql')) {
      return "$input.command < ${args[0]} > '$input.produces'"
    }
    return "$input.command ${args.join(' ')} '$input.produces'"
  }

  // String scripts(List<UnmodifiableMap> inputs) {
  //   List<String> scripts = inputs.inject([]) { acc, inp -> acc << Input.script(inp) }
  //   return scripts.join('\n')
  // }

  static Boolean is_excluded(Map input, String filename) {
    return input.exclude.inject(false) { agg, pattern -> agg || (filename =~ pattern) }
  }
}

