#!/bin/bash

if [ -z "$VIRTUAL_ENV" ]; then
  echo "Activate venv."
  cd ..
  source venv/bin/activate
  cd lab3
fi

cd moon_gen
if [ ! -d "build" ]; then
  cmake -B build
fi

cmake --build build 
build/moon

cd ..
if [ ! -d "build" ]; then
  cmake -B build
fi
cmake --build build
ffc -l dolfin Poisson.ufl
dolfin-convert moon.msh moon.xml
build/demo_poisson