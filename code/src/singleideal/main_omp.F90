program streams_singleideal_omp

use streams_equation_singleideal_omp_object
use streams_parameters
use mpi

implicit none
type(equation_singleideal_omp_object) :: singleideal_omp
integer :: mpi_err

call mpi_initialize()

call singleideal_omp%run(filename='singleideal.ini')

call MPI_Finalize(mpi_err)

endprogram streams_singleideal_omp

