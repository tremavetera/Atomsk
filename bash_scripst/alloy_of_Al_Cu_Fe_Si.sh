#!/bin/bash

atomsk --create fcc 4.046 Al orient [100] [010] [001] \
-duplicate 2 2 2 \
-select random 25% Al \
-substitute Al Fe \
-select random 33.33% Al \
-substitute Al Si \
-select random 50% Al \
-substitute Al Cu \
AL25Fe25Si25Cu25.cif
