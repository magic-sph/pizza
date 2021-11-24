file(GLOB srcfiles "*.f90")
list(REMOVE_ITEM srcfiles "${SRC}/pizza.f90")

set_source_files_properties(mpimod.f90 PROPERTIES COMPILE_FLAGS ${FFLAG_NO_STD})

add_executable(pizza.exe pizza.f90 ${srcfiles})
set_property(TARGET ${EXEC} PROPERTY LINKER_LANGUAGE Fortran)

# Link MPI
if(USE_MPI MATCHES yes)
   if(MPI_Fortran_FOUND AND NOT MPI_WRAPPER_FOUND)
      target_link_libraries(${EXEC} ${MPI_Fortran_LIBRARIES})
   endif()
endif()

if(USE_LAPACKLIB MATCHES MKL)
   add_library(lapack95 ${MKLROOT}/include/lapack.f90)
   add_library(blas95 ${MKLROOT}/include/blas.f90)
   target_link_libraries(lapack95 ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
   target_link_libraries(blas95 ${MKLCORE} ${MKLBLAS} ${MKLSEQ} ${MKL64})
   target_link_libraries(${EXEC} lapack95 blas95)
endif()

# Link LAPACK
if(USE_LAPACKLIB MATCHES LAPACK)
   target_link_libraries(${EXEC} ${LAPACK_LIBRARIES})
elseif(USE_LAPACKLIB MATCHES LIBSCI)
   target_link_libraries (${EXEC} ${SCI})
endif()

target_link_libraries(${EXEC} ${FFTW3})
