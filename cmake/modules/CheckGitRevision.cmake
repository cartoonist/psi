# Check git revision.
#
# This set the following variables if git information is available:
#   - GIT_REVISION
#   - GIT_DESC
#   - GIT_BRANCH
#   - GIT_COMMIT_ID
#   - GIT_COMMIT_DATE

# Get the git revision
execute_process(COMMAND git describe --always
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_DESC
  ERROR_QUIET)
string(STRIP "${GIT_DESC}" GIT_DESC)

# Check whether we got any revision
if (NOT "${GIT_DESC}" STREQUAL "")
  execute_process(COMMAND git name-rev --name-only HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH)
  execute_process(COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_ID)
  execute_process(COMMAND git log -1 --format=%cd --date=format:%h\ %d\ %Y
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DATE)
  string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
  string(REGEX MATCH "[^/]*$" GIT_BRANCH ${GIT_BRANCH})
  string(STRIP "${GIT_COMMIT_ID}" GIT_COMMIT_ID)
  string(STRIP "${GIT_COMMIT_DATE}" GIT_COMMIT_DATE)
  if ("${GIT_DESC}" STREQUAL "v${PROJECT_VERSION}"
      OR "${GIT_DESC}" STREQUAL "${PROJECT_VERSION}")
    set(GIT_REVISION "${GIT_DESC}")
  else()
    set(GIT_REVISION "v${PROJECT_VERSION}-${GIT_BRANCH}-${GIT_COMMIT_ID}")
  endif()
endif()
