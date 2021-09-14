rm -r make
mkdir make
cd make
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
#echo 'SET (CMAKE_C_FLAGS_INIT                "-Wall -Wextra -Werror -std=c99")
echo 'SET (CMAKE_C_FLAGS_INIT                "-Wall -std=c99")
SET (CMAKE_C_FLAGS_DEBUG_INIT          "-g")
SET (CMAKE_C_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
SET (CMAKE_C_FLAGS_RELEASE_INIT        "-O3 -DNDEBUG")
SET (CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

SET (CMAKE_CXX_FLAGS_INIT                "-Wall")
SET (CMAKE_CXX_FLAGS_DEBUG_INIT          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE_INIT        "-O3 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")' >> ~/ClangOverrides.txt

cmake -DCMAKE_USER_MAKE_RULES_OVERRIDE=~/ClangOverrides.txt -D_CMAKE_TOOLCHAIN_PREFIX=llvm- ../src
make -j4
make install
