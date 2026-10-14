cmake_minimum_required(VERSION 3.15)   # 4.3 yerine daha yaygın bir sürüm
project(simulasyon)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(simulasyon ins.cpp)     # SADECE bu dosyayı ekleyin