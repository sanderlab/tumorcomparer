#!/bin/bash

echo "Pulling new version"
sudo docker pull cannin/tumorcomparer

echo "Stopping new container"
sudo docker stop tumorcomparer
sudo docker rm tumorcomparer

echo "Starting new container"
sudo docker run --restart always --name tumorcomparer -d -p 3843:3838 -t cannin/tumorcomparer shiny-server
