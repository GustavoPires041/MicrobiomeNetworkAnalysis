#!/bin/bash

if [ ! -d "reports/clusters/output" ]; then
    echo "A pasta 'output' não existe. Criando..."
    mkdir -p reports/clusters/output
else
    echo "A pasta 'output' já existe."
fi


for file in reports/clusters/input/*; do
    echo "$file"
    if [ -f "$file" ]; then

        filename=$(basename -- "$file")
        filename_no_ext="${filename%.*}"


        manta -i "$file" -o "reports/clusters/output/${filename_no_ext}_cluster" -b -seed 1412
    fi
done
