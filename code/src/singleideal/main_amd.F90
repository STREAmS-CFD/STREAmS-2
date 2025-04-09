program streams_singleideal_amd

use streams_equation_singleideal_amd_object
use streams_parameters
use mpi

implicit none
type(equation_singleideal_amd_object) :: singleideal_amd
integer :: mpi_err

call mpi_initialize()

call singleideal_amd%run(filename='singleideal.ini')

call MPI_Finalize(mpi_err)

endprogram streams_singleideal_amd

