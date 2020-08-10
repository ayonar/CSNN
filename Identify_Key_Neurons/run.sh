#!/bin/bash

for f in *.ipynb
do
    jupyter nbconvert \
      --ExecutePreprocessor.allow_errors=True \
      --ExecutePreprocessor.timeout=-1 \
      --FilesWriter.build_directory=../results \
      --execute "$f" &

wait
done

mkdir -p ../results/pkl && mv *.pkl ../results/pkl