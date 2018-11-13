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
    def data = new LinkedHashMap(Parse.defaults(db_name));
    if (!spec) {
      return data;
    }
    spec.directives.each { entry ->
      data.directives[entry.key] = entry.value;
    }
    data += spec.subMap('command', 'produces');
    return data;
  }
}
