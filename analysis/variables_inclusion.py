from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists


def generate_inclusion_variables(index_date):
  inclusion_variables = dict(
    
    registered = patients.registered_as_of(
        index_date,
    ), 

    has_died = patients.died_from_any_cause(
      on_or_before=index_date,
      returning="binary_flag",
    ),
          
  )
  return inclusion_variables

