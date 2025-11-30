# Базовое использование
mpirun -np 4 ./petsc_solver -n 1000

# С различными предобуславливателями
mpirun -np 4 ./petsc_solver -n 2000 -pc_type ilu

# Запуск тестов
make test

# Запуск примеров
mpirun -np 4 examples/poisson2d -nx 100 -ny 100