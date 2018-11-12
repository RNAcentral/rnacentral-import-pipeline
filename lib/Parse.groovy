class Parse {
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

  static Map build(String db_name, Map spec) {
    def defaults = Parse.defaults(db_name);
    return defaults + (spec ? spec : [:]);
  }
}
