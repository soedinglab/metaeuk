include_directories(commons)
add_subdirectory(assembler)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(plass
        ${commons_source_files}
        ${assembler_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        plass.cpp
        )
mmseqs_setup_derived_target(plass)
target_link_libraries(plass version kerasify)

add_dependencies(plass local-generated)

install(TARGETS plass DESTINATION bin)
