program streams_singleideal_cpu

use streams_equation_singleideal_cpu_object
use streams_parameters
use mpi

implicit none
type(equation_singleideal_cpu_object) :: singleideal_cpu
integer :: mpi_err

call mpi_initialize()

call singleideal_cpu%run(filename='singleideal.ini')

call MPI_Finalize(mpi_err)

endprogram streams_singleideal_cpu

