program streams_singleideal_gpu
!< STREAmS, STREAmS for Navier-Stokes equation, GPU backend.

use streams_equation_singleideal_gpu_object
use streams_parameters
use mpi

implicit none
type(equation_singleideal_gpu_object) :: singleideal_gpu !< Navier-Stokes equations system.
integer :: mpi_err

call mpi_initialize()

call singleideal_gpu%run(filename='singleideal.ini')

call MPI_Finalize(mpi_err)

endprogram streams_singleideal_gpu

