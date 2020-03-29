$build_folder = "./vsproject"
if(!(Test-Path $build_folder)) {
    md $build_folder
}
cd $build_folder
cmake ../ -G "Visual Studio 16 2019" -A x64 -DCMAKE_CONFIGURATION_TYPES="Release;Debug"
cd ..