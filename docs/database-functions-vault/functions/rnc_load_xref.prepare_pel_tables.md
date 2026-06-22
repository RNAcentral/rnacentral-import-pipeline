# rnc_load_xref.prepare_pel_tables

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`prepare_pel_tables__66815.sql`](../../../database_functions/rnc_load_xref/rnc_load_xref/prepare_pel_tables__66815.sql)

```sql
prepare_pel_tables() RETURNS void
```

(Re)creates the build tables `xref_pel_deleted` and `xref_pel_not_deleted` as `LIKE xref`
(defaults/constraints/comments/generated, **without** indexes/FKs), and drops any leftover
indexes and FK constraints so they load fast. Step 1 of [[rnc_load_xref.load_xref]].

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `xref_pel_deleted`, `xref_pel_not_deleted`, `xref`
