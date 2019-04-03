file(GLOB srcfiles "*.f90")
list(REMOVE_ITEM srcfiles "${SRC}/pizza.f90")

set_source_files_properties(mpimod.f90 PROPERTIES COMPILE_FLAGS ${FFLAG_NO_STD})

if(USE_SHTNS MATCHES yes AND NOT SHTNS MATCHES SHTNS-NOTFOUND)
else()
   list(REMOVE_ITEM srcfiles "${SRC}/shtns.f90")
endif()


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
endif()

target_link_libraries(${EXEC} ${FFTW3})

if(USE_SHTNS MATCHES yes AND NOT SHTNS MATCHES SHTNS-NOTFOUND)
   # check shtns backend for ffts (mkl or fftw)
   FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90"
   "
   program TestFFTW
    call shtns_verbose(0)
   end program TestFFTW
   ")

   try_run(SHTNS_RUN SHTNS_LINK ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      LINK_LIBRARIES ${SHTNS} ${MKLCORE} ${MKLSEQ} ${MKL64}
      COMPILE_OUTPUT_VARIABLE OUTPUT)

   if(SHTNS_LINK)
      set(SHTNS_DEPS ${MKLCORE} ${MKLSEQ} ${MKL64})
   endif()

   if(NOT SHTNS_LINK)
      set(SAFE_CMAKE_LIBS "${LIBS}")
      file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90"
      "
      program TestFFTW
      call shtns_verbose(0)
      end program TestFFTW
      ")
      try_run(SHTNS_RUN SHTNS_LINK ${CMAKE_BINARY_DIR}
         ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90
         COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
         LINK_LIBRARIES ${SHTNS} ${FFTW3_OMP} ${FFTW3}
         COMPILE_OUTPUT_VARIABLE OUTPUT)
      if(SHTNS_LINK)
         set(SHTNS_DEPS ${FFTW3_OMP} ${FFTW3})
      endif(SHTNS_LINK)
   endif(NOT SHTNS_LINK)

   set(LIBS "${LIBS} ${SHTNS_DEPS}")
   target_link_libraries(${EXEC} ${SHTNS})
   target_link_libraries(${EXEC} ${SHTNS_DEPS})

endif()
