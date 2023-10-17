# How to add something to search export



### Rust

In `utils/search-export/src/main.rs` Add

1. A new member of the `Groupable` enumeration. It should have a good name that tells you what it is
2. A new cli argument. The new argument should be placed before the output argument, to facilitate future expansion easily.
3. In the `main` function, add a new branch to the match to handle the new argument. It should dispatch to the new function you will write to handle the new search facet
4. In the call to `sequences::writers::write_merge`, add the new search facet. It can go anywhere.

In a new file under the correct subdirectory (`genes` or `sequences`) add:

1. A struct definition with the right members and a public `id` member of type `usize`. Derive from serde `Deserialize` and `Serialize`, `Clone`, `Debug`, `PartialEq` and `Eq`.
2. An implementation for `grouper::HasIndex` for your struct.
3. A `group` function.
4. Getter functions for the private members of the struct.

See the other export code for examples and use them like a template.

In e.g. `utils/search-export/src/sequences/file_joiner.rs`

1. Add your new facet to the list of things in `use::super`
2. Add your new term to the enumeration of `FileTypes`
3. Add your new term to the FileJoiner struct - see the others for an example, but it has to be a StreamDeserializer
4. Create the iterator for the new term in `pub fn build`, and include it in the returned value in the right place
5. Add the returned iterator to the implementation of `next for FileJoiner`. You will need to add some things to the `match` statement as well, you should be able to figure it out from the others

In e.g. `utils/search-export/src/sequences/raw.rs`

1. Add your new facet to the `use crate::sequences::{` list
2. Add your new facet to the `pub struct Raw`. The exact type will depend on what type of grouping you used, e.g. Only one -> bare type, Multiple -> Vec\<type\> and optional -> Some(type)
3. A public getter function to get a reference to the new facet data from `Raw`

In e.g. `utils/search-export/src/sequences/normalized.rs`

1. Include your new type in the `use crate::sequences` list
2. Add your new type to the `pub struct Normalized`, It will be of whatever type you used in Raw, e.g. bare type, Option(type) or Vec\<type\>
3. Add the new type to the construction wrapped in `Ok(` at the end opf the `pub fn new` implementation for `Normalized`

Don't forget to add the new module in the `mod.rs` file in whichever subdirectory you used

### Python

In `rnacentral_pipeline/rnacentral/search_export/data.py` add:

1. A getter function that converts the value exported from the database into the right type for search export (e.g. into a bool or something)
2. A new field in the call to `entry` to create `builder`. Your new search field will go in the call to `section` in the list of `"additional_fields"`. **DO NOT PUT IT AT THE END OF THE LIST!** The last thing in the list should be the heirachical type for the SO type tree. Putting a standard field after that will make modifying the xml really difficult if it is necessary.
3.

### Nextflow
The correct bit of the workflow to insert into depends on the thing being exported for.

In e.g. `sequences.nf` Add

1. A process to copy the search data out of the database and apply the rust `search-export group <new facet>` code to it.
2. A new channel from the relevant SQL file
3. A call to the process inside the `build_metadata` call
4. The correct argument in the right place in the `build_metadata` definition and script. You'll need a new `path` input to the nextflow process, and to put the argument in the right place for the cli call to the rust code.

### SQL
Only one thing to add, and it completely depends on the facet you're adding. You need

1. A SQL file that extracts: 'search_export_urs.id`, 'search_export_urs.urs_taxid` and whatever data goes into the search facet. This will almost certainly require a join from `search_export_urs` to some other table containing the data you would like to be searchable.
