#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { select } from './workflows/databases/select.nf'


workflow {
  emit: done

  select() | set { done }

}
