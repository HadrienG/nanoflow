#!/usr/bin/env bash

cd $HOME
curl -fsSL get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
