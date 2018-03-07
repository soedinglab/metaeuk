include_directories(commons)
add_subdirectory(exonpredictor)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(metaeuk
        ${commons_source_files}
        ${exonpredictor_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        metaeuk.cpp
        )
mmseqs_setup_derived_target(metaeuk)
target_link_libraries(metaeuk version)

add_dependencies(metaeuk local-generated)

install(TARGETS metaeuk DESTINATION bin)
