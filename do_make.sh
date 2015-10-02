#!/bin/bash

cmake -GNinja  -DDEAL_II_DIR="/u/qiaolei/devel/deal.II/testWithoutInstall/build" .
ninja run
