# CMake generated Testfile for 
# Source directory: /home/bk/prog/eng/thermalFVM
# Build directory: /home/bk/prog/eng/thermalFVM
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(thermal "/home/bk/prog/eng/thermalFVM/UnitTests/thermal_test")
ADD_TEST(radiation "/home/bk/prog/eng/thermalFVM/UnitTests/radiation_test")
ADD_TEST(HTCDittusBoelter "/home/bk/prog/eng/thermalFVM/UnitTests/HTCDittusBoelter_test")
ADD_TEST(HTC "/home/bk/prog/eng/thermalFVM/UnitTests/HTC_test")
ADD_TEST(temperature "/home/bk/prog/eng/thermalFVM/UnitTests/temperature_test")
ADD_TEST(heatFlux "/home/bk/prog/eng/thermalFVM/UnitTests/heatFlux_test")
