nextflow.enable.dsl=2

process fetch_training {
  input:
  path(query)

  output:
  path('data.csv'), emit: train
  path('test.csv', emit: test

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > raw
  rnac genes normalize-training --test-split $params.genes.test_split raw train.csv test.csv
  """
}

process train {
  input:
  path(data)

  output:
  path('genes-model.joblib')

  """
  rnac genes train $data genes-model.joblib
  """
}

process validate {
  input:
  tuple path(model), path(test_data)

  output:
  tuple path(model), path('genes-validation.log')

  """
  rnac genes validate-model $model $test_data genes-validation.log
  """
}

process save_work {
  input:
  tuple path(model), path(log)

  output:
  path(model)

  """
  cp $model $params.genes.model
  cp $log $params.genes.validation
  """
}


workflow train {
  emit:
    model
  main:
    Channel.ofPath('files/genes/training-data.sql') | fetch_training

    fetch_training.out.train | train | set { trained }

    trained.combine(fetch_training.out.test) \
    | validate \
    | save_work \
    | set { model }
}

workflow {
  train()
}
