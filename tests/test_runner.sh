#!/bin/bash

set -e

echo "Running PETSc GMRES solver tests..."

# Компиляция тестов
make test

# Запуск тестов
mpirun -np 2 ./test_solver -ksp_monitor -ksp_converged_reason

echo "All tests completed successfully!"