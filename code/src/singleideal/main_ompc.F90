program streams_singleideal_ompc

use streams_equation_singleideal_ompc_object
use streams_parameters
use mpi

implicit none
type(equation_singleideal_ompc_object) :: singleideal_ompc
integer :: mpi_err

call mpi_initialize()

call singleideal_ompc%run(filename='singleideal.ini')

call MPI_Finalize(mpi_err)

endprogram streams_singleideal_ompc

