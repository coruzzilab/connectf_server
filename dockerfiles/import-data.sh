#!/usr/bin/env bash

: '
This script runs the import commands specified in the readme file 
for the specified IMPORT directory: arabidopsis, maize or rice:

  python manage.py import_annotation -i <annotation_file>
  python manage.py import_data <data_dir> <metadata_dir>
  python manage.py import_edges <edges_file>
'

if [ -n "$IMPORT" ]; then
  ANNOTATION_FILE=$(find ./connectf_data_release_v1/$IMPORT -name "annotation*.csv.gz")

  # Run migrations
  python manage.py migrate

  # Look for annotation file, if found, import it
  if [ -f "$ANNOTATION_FILE" ]; then
    echo "Importing annotation file: $ANNOTATION_FILE"
    python manage.py import_annotation -i $ANNOTATION_FILE
  else
    echo "Annotation file not found: $ANNOTATION_FILE"
  fi

  # Import data
  echo "Importing data from ./connectf_data_release_v1/$IMPORT"
  python manage.py import_data ./connectf_data_release_v1/$IMPORT/data ./connectf_data_release_v1/$IMPORT/metadata

  # for each file either .csv or .csv.gz found in ./connectf_data_release_v1/$IMPORT/additional_edges, import it
  for file in ./connectf_data_release_v1/$IMPORT/additional_edges/*; do
    if [ -f "$file" ]; then
      python manage.py import_edges $file
    fi
  done
fi