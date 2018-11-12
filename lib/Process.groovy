class Process {
  static Map defaults(db_name) {
    return [
      command: "rnac external $db_name",
      produces: '*.csv',
      directives: [
        memory: 2.GB,
        tag: db_name,
        group: "/rnc/process/$db_name"
      ]
    ];
  }

  static Map build(db_name, spec) {
    def defaults = Process.defaults(db_name);
    return defaults + spec;
  }

}
