# Practica-MPI
Desarrollo de una versión paralelizada con [MPI-2.0](https://www.open-mpi.org/software/ompi/v2.0/) del programa propuesto.

MPI (*Message Passing Interface*) es un estándar de paso de mensajes estandarizado y portátil diseñado para funcionar en arquitecturas de computación paralela. El estándar MPI define la sintaxis y la semántica de las rutinas de la biblioteca que son útiles para una amplia gama de usuarios que escriben programas portátiles de paso de mensajes en C, C++ y Fortran. Su principal característica es que no precisa de memoria compartida, por lo que es muy importante en la programación de sistemas distribuidos.

### Compilar un programa MPI
```c
mpicc programa.c -o programa
```

### Ejecutar un programa MPI
```c
mpiexec -n <numprocs> ./programa
```
