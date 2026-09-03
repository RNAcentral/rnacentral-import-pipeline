set search_path = rnacen, public;
set client_min_messages = warning;

-- sequences referenced in rna (fk4 target)
insert into rnacen.rna (urs, len) values
  ('UPI_A',100),('UPI_A2',100),('UPI_B',100),('UPI_C',100),
  ('UPI_D',100),('UPI_E',100),('UPI_F',100);

-- Initial xref state = result of previous release (id 1). created=last=1.
insert into rnacen.xref (ac, dbid, version, version_i, urs, created, last, deleted, taxid) values
  ('ACC1', 1, 1, 1, 'UPI_A', 1, 1, 'N', NULL),   -- unchanged (also gains a new variant)
  ('ACC2', 1, 1, 1, 'UPI_B', 1, 1, 'N', NULL),   -- version will change (case 2)
  ('ACC3', 1, 1, 1, 'UPI_C', 1, 1, 'N', NULL),   -- dropped from load (case 5)
  ('ACC5', 1, 1, 1, 'UPI_E', 1, 1, 'Y', NULL),   -- previously deleted -> reactivated (case D)
  ('ACC6', 1, 1, 1, 'UPI_F', 1, 1, 'Y', NULL);   -- stays deleted (carried forward)

-- Incoming load for release 2. comparable_prot_upi already resolved to the urs.
insert into rnacen.load_retro_tmp (in_dbid, in_load_release, in_ac, in_version, in_taxid, comparable_prot_upi) values
  (1, 2, 'ACC1', 1, 9606, 'UPI_A'),    -- unchanged   -> refresh
  (1, 2, 'ACC1', 1, 9606, 'UPI_A2'),   -- new variant -> Gap A (version_i = 2)
  (1, 2, 'ACC2', 2, 9606, 'UPI_B'),    -- version bump-> retire old + new gen (version_i 1)
  (1, 2, 'ACC4', 1, 9606, 'UPI_D'),    -- brand new   -> case D (version_i 1)
  (1, 2, 'ACC5', 1, 9606, 'UPI_E');    -- reactivate  -> case D (version_i 1)
