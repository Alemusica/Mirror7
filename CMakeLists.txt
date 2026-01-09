cmake_minimum_required(VERSION 3.22)
project(Mirror7 VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

# Prefer reproducible macOS builds out of the box
if(APPLE)
    if(NOT CMAKE_OSX_ARCHITECTURES)
        set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64" CACHE STRING "" FORCE)
    endif()
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "" FORCE)
    endif()
endif()

# Root of the aureonoise mono-repo (used for shared DSP headers and spatial assets)
set(AUREONOISE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../aureonoise_tilde" CACHE PATH "Path to the aureonoise mono-repo")
set(AUREONOISE_CORE_DIR "${AUREONOISE_ROOT}/source/aureo_core" CACHE PATH "Path to aureo_core headers")
set(AUREONOISE_LEGACY_DIR "${AUREONOISE_ROOT}/source/aureonoise3beta/legacy" CACHE PATH "Path to beta7_tools legacy sources")

if(NOT EXISTS "${AUREONOISE_CORE_DIR}/aureo_noise.hpp")
    message(FATAL_ERROR "aureo_core headers not found. Set AUREONOISE_ROOT to the repo containing source/aureo_core.")
endif()

set(_MIRROR7_REQUIRED_PHI "${AUREONOISE_LEGACY_DIR}/beta7_tools/include/aureonoise/spatial/phi_model.hpp")
if(NOT EXISTS "${_MIRROR7_REQUIRED_PHI}")
    set(_MIRROR7_LEGACY_FALLBACKS
        "${AUREONOISE_ROOT}/source/vellutoblu~"
        "${AUREONOISE_ROOT}/vellutoblu~"
    )
    foreach(_candidate IN LISTS _MIRROR7_LEGACY_FALLBACKS)
        if(EXISTS "${_candidate}/beta7_tools/include/aureonoise/spatial/phi_model.hpp")
            set(AUREONOISE_LEGACY_DIR "${_candidate}" CACHE PATH "Path to beta7_tools legacy sources" FORCE)
            set(_MIRROR7_REQUIRED_PHI "${AUREONOISE_LEGACY_DIR}/beta7_tools/include/aureonoise/spatial/phi_model.hpp")
            message(STATUS "Mirror7: using legacy modules from ${AUREONOISE_LEGACY_DIR}")
            break()
        endif()
    endforeach()
endif()

if(NOT EXISTS "${_MIRROR7_REQUIRED_PHI}")
    message(FATAL_ERROR "beta7_tools spatial headers not found under AUREONOISE_LEGACY_DIR (${AUREONOISE_LEGACY_DIR}). Set AUREONOISE_LEGACY_DIR manually.")
endif()
unset(_MIRROR7_LEGACY_FALLBACKS)
unset(_MIRROR7_REQUIRED_PHI)

find_package(JUCE CONFIG REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(MIRROR7_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/source)
set(MIRROR7_RESOURCES_DIR ${CMAKE_CURRENT_LIST_DIR}/resources)
set(AUREONOISE_SPATIAL_ASSETS_DIR "${MIRROR7_RESOURCES_DIR}/spatial_assets" CACHE PATH "Location of spatial assets")

add_library(mirror7_engine STATIC
    ${MIRROR7_SOURCE_DIR}/engine/Mirror7Engine.cpp
    ${MIRROR7_SOURCE_DIR}/engine/DialogueSystem.cpp
    ${MIRROR7_SOURCE_DIR}/engine/NoiseController.cpp
    ${MIRROR7_SOURCE_DIR}/engine/Spatializer.cpp
    ${AUREONOISE_LEGACY_DIR}/beta7_tools/src/spatial/phi_model.cpp
    ${AUREONOISE_LEGACY_DIR}/beta7_tools/src/spatial/cefg.cpp
)

target_include_directories(mirror7_engine
    PUBLIC
        ${MIRROR7_SOURCE_DIR}
        ${AUREONOISE_CORE_DIR}
        ${AUREONOISE_ROOT}/source
        ${AUREONOISE_LEGACY_DIR}
        ${AUREONOISE_LEGACY_DIR}/core
        ${AUREONOISE_LEGACY_DIR}/beta7_tools/include
)

target_link_libraries(mirror7_engine
    PUBLIC
        juce::juce_dsp
)

target_compile_definitions(mirror7_engine
    PUBLIC
        AUREONOISE_SPATIAL_ASSETS_DIR="${AUREONOISE_SPATIAL_ASSETS_DIR}"
)

juce_add_plugin(Mirror7
    PRODUCT_NAME "Mirror7"
    COMPANY_NAME "Aureonoise"
    PLUGIN_MANUFACTURER_CODE AuNo
    PLUGIN_CODE Mir7
    FORMATS VST3 AU Standalone
    IS_SYNTH TRUE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    COPY_PLUGIN_AFTER_BUILD TRUE
    VST3_CAN_REPLACE_VST2 FALSE
)

target_sources(Mirror7
    PRIVATE
        ${MIRROR7_SOURCE_DIR}/plugin/PluginProcessor.cpp
        ${MIRROR7_SOURCE_DIR}/plugin/PluginProcessor.h
        ${MIRROR7_SOURCE_DIR}/gui/PluginEditor.cpp
        ${MIRROR7_SOURCE_DIR}/gui/PluginEditor.h
        ${MIRROR7_SOURCE_DIR}/config/ParamIDs.h
)

target_include_directories(Mirror7
    PRIVATE
        ${MIRROR7_SOURCE_DIR}
)

target_compile_definitions(Mirror7
    PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
)

target_link_libraries(Mirror7
    PRIVATE
        mirror7_engine
        juce::juce_audio_utils
        juce::juce_audio_processors
        juce::juce_gui_extra
)

target_compile_definitions(Mirror7 PRIVATE
    AUREONOISE_SPATIAL_ASSETS_DIR="${AUREONOISE_SPATIAL_ASSETS_DIR}"
)

target_include_directories(Mirror7 PRIVATE ${AUREONOISE_CORE_DIR})

# Post-build moduleinfo fixer (VST3 manifest cleanup)
set(MIRROR7_ARTEFACT_ROOT "${CMAKE_CURRENT_BINARY_DIR}/Mirror7_artefacts")
set(MIRROR7_FIXER_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/scripts/post_build_fix_moduleinfo.py")

add_custom_command(TARGET Mirror7 POST_BUILD
    COMMAND ${Python3_EXECUTABLE} "${MIRROR7_FIXER_SCRIPT}"
            --artefact-root "${MIRROR7_ARTEFACT_ROOT}"
            --vst3-name "Mirror7"
)

if(TARGET Mirror7_VST3)
    add_custom_command(TARGET Mirror7_VST3 POST_BUILD
        COMMAND ${Python3_EXECUTABLE} "${MIRROR7_FIXER_SCRIPT}"
                --artefact-root "${MIRROR7_ARTEFACT_ROOT}"
                --vst3-name "Mirror7"
    )
endif()

if(BUILD_TESTING)
    add_executable(mirror7_engine_smoke
        tests/EngineSmokeTest.cpp
    )
    target_link_libraries(mirror7_engine_smoke
        PRIVATE
            mirror7_engine
    )
    add_test(NAME mirror7_engine_smoke COMMAND mirror7_engine_smoke)

    add_executable(mirror7_dialogue_test
        tests/DialogueTest.cpp
    )
    target_link_libraries(mirror7_dialogue_test
        PRIVATE
            mirror7_engine
    )
    add_test(NAME mirror7_dialogue_test COMMAND mirror7_dialogue_test)

    add_executable(mirror7_spatial_test
        tests/SpatializerTest.cpp
    )
    target_link_libraries(mirror7_spatial_test
        PRIVATE
            mirror7_engine
    )
    add_test(NAME mirror7_spatial_test COMMAND mirror7_spatial_test)
endif()
